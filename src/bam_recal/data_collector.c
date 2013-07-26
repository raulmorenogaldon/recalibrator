#include "data_collector.h"
#include "string.h"

static char ult_seq[100];
static int l_ult_seq = 0;
static int pos_ult_seq;

/**
 * Get recalibration data from BAM path.
 */
ERROR_CODE
recal_get_data_from_file(const char *bam_path, const char *ref_name, const char *ref_path, recal_info_t *out_info)
{
	genome_t* ref;
	bam_file_t *bam_f;

	//Open bam
	printf("Opening BAM from \"%s\" ...\n", bam_path);
	bam_f = bam_fopen(bam_path);
	printf("BAM opened!...\n");

	//Open reference genome
	printf("Opening reference genome from \"%s%s\" ...\n", ref_path, ref_name);
	ref = genome_new(ref_name, ref_path);
	printf("Reference opened!...\n");

	//Fill data
	#ifdef D_TIME_DEBUG
		time_init_slot(D_SLOT_GET_DATA_BAM, clock(), TIME_GLOBAL_STATS);
	#endif
	recal_get_data_from_bam(bam_f, ref, out_info);
	#ifdef D_TIME_DEBUG
		time_set_slot(D_SLOT_GET_DATA_BAM, clock(), TIME_GLOBAL_STATS);
	#endif

	//Memory free
	printf("\nClosing BAM file...\n");
	bam_fclose(bam_f);
	genome_free(ref);
	printf("BAM closed.\n");

	return NO_ERROR;
}

/**
 * Get recalibration data from BAM file.
 */
ERROR_CODE
recal_get_data_from_bam(const bam_file_t *bam, const genome_t* ref, recal_info_t* output_data)
{
	bam_batch_t* batch;
	uint8_t *last_seq;
	uint32_t l_last_seq;
	uint32_t pos_last_seq;
	bam1_t *last_alig;

	int count = 0;

	#ifdef USE_BATCH_POOL
	bam_pool_t* pool;
	#endif

	int zero = 0;

	unmapped = 0;
	duplicated = 0;
	#ifdef NOT_MAPPING_QUAL_ZERO
	mapzero = 0;
	#endif

	//DEBUG
	#ifdef D_SAMPLES_OUTPUT
	FILE *fp;
	out_path = "cadenas.out";
	fp = fopen(out_path, "w");
	fclose(fp);
	#endif

	printf("\n----------------\nProcessing \"%s\" file...\n----------------", bam->filename);

	//Allocate memory for batchs
	batch = bam_batch_new(MAX_BATCH_SIZE, MULTIPLE_CHROM_BATCH);

	//Read batch
	#ifdef D_TIME_DEBUG
		time_init_slot(D_SLOT_READ_BATCH, clock(), TIME_GLOBAL_STATS);
	#endif
	bam_fread_max_size(batch, MAX_BATCH_SIZE, 0, bam);
	#ifdef D_TIME_DEBUG
		time_set_slot(D_SLOT_READ_BATCH, clock(), TIME_GLOBAL_STATS);
	#endif

	printf("\nNum alignments in batchs: %d\n----------------\n", batch->num_alignments);

	#ifdef USE_BATCH_POOL
	printf("Using batch pool!\n");
	pool = pool_new(3, MAX_BATCH_SIZE);
	#endif

	last_seq = (uint8_t *)malloc(sizeof(uint8_t));	// Avoid comprobation in bucle and always use free
	while(batch->num_alignments != 0
		#ifdef D_MAX_READS
			&& count < D_MAX_READS
		#endif
		)
	{
		#ifdef USE_BATCH_POOL
		while(!pool_has_free_batchs(pool));
		{
			#ifdef D_TIME_DEBUG
				time_init_slot(D_SLOT_MEMCOPY_BATCH, clock(), time_global_stats);
			#endif
			pool_add_batch(batch, pool);
			bam_batch_free(batch, 1);

			//Get batch in pool
			pool_next_batch(pool);
			batch = pool_get_batch(pool);
			#ifdef D_TIME_DEBUG
				time_set_slot(D_SLOT_MEMCOPY_BATCH, clock(), time_global_stats);
			#endif
		}
		#endif

		//Process batch
		#ifdef D_TIME_DEBUG
			time_init_slot(D_SLOT_GET_DATA_BATCH, clock(), TIME_GLOBAL_STATS);
		#endif
		recal_get_data_from_bam_batch(batch, ref, output_data);
		#ifdef D_TIME_DEBUG
			time_set_slot(D_SLOT_GET_DATA_BATCH, clock(), TIME_GLOBAL_STATS);
		#endif

		//Update read counter
		count += batch->num_alignments;

		//Show total progress
		printf("\n Total alignments readed: %d", count);

		//Get last alignment
		free(last_seq);
		last_alig = batch->alignments_p[batch->num_alignments-1];
		last_seq = new_sequence_from_bam(last_alig);
		l_last_seq = last_alig->core.l_qseq;
		pos_last_seq = last_alig->core.pos;

		//Free memory and take a new batch
		#ifndef USE_BATCH_POOL
		bam_batch_free(batch, 1);
		#endif
		batch = bam_batch_new(MAX_BATCH_SIZE, MULTIPLE_CHROM_BATCH);

		//Read batch
		#ifdef D_TIME_DEBUG
			time_init_slot(D_SLOT_READ_BATCH, clock(), TIME_GLOBAL_STATS);
		#endif
		//bam_fread_max_size(batch, MAX_BATCH_SIZE, 1, bam);
		bam_fread_max_size_no_duplicates(batch, MAX_BATCH_SIZE, 0, bam, last_seq, &l_last_seq, &pos_last_seq);
		#ifdef D_TIME_DEBUG
			time_set_slot(D_SLOT_READ_BATCH, clock(), TIME_GLOBAL_STATS);
		#endif
	}

	printf("\n----------------\n%d alignments readed.", count);
	printf("\n%d alignments processed.", count - unmapped - mapzero - duplicated);
	printf("\n%d alignments duplicated.", duplicated);
	#ifdef NOT_MAPPING_QUAL_ZERO
	printf("\n%d alignments with map quality zero.", mapzero);
	#endif
	printf("\n%d alignments unmapped.", unmapped);

	#ifdef USE_BATCH_POOL
		pool_free(pool);
	#endif

	//Free batch
	bam_batch_free(batch, 1);
	free(last_seq);

	return NO_ERROR;
}

/**
 * Get recalibration data from BAM batch of alignments.
 */
ERROR_CODE
recal_get_data_from_bam_batch(const bam_batch_t* batch, const genome_t* ref, recal_info_t* output_data)
{
	int i;

	//Process all alignments of the batchs
	for(i = 0; i < batch->num_alignments; i++)
	{
		//Process every alignment
		#ifdef NOT_MAPPING_QUAL_ZERO
		if(batch->alignments_p[i]->core.qual != 0)
		#endif
		{
		#ifdef D_TIME_DEBUG
			time_init_slot(D_SLOT_GET_DATA_ALIG, clock(), TIME_GLOBAL_STATS);
		#endif
		recal_get_data_from_bam_alignment(batch->alignments_p[i], ref, output_data);
		#ifdef D_TIME_DEBUG
			time_set_slot(D_SLOT_GET_DATA_ALIG, clock(), TIME_GLOBAL_STATS);
		#endif
		}
		#ifdef NOT_MAPPING_QUAL_ZERO
		else
		{
			mapzero++;
		}
		#endif
	}

	return NO_ERROR;
}

/**
 * Get recalibration data from alignment.
 */
ERROR_CODE
recal_get_data_from_bam_alignment(const bam1_t* alig, const genome_t* ref, recal_info_t* output_data)
{
	alignment_t *aux_alig;
	uint8_t *quals;
	uint8_t *bam_seq;
	uint8_t *ref_seq;
	uint8_t aux_comp[16];
	uint32_t init_pos, end_pos;
	uint32_t cycles;
	uint8_t *comp_res;
	uint8_t *dinucs;
	uint32_t flag;

	#ifdef __SSE2__
	__m128i v_ref, v_seq, v_comp;
	#endif

	unsigned int i, j;

	//Unmapped readings
	if(alig->core.tid < 0)
	{
		unmapped++;
		return NO_ERROR;
	}

	//aux_alig = alignment_new_by_bam(alig, 1);

	//Get sequence
	bam_seq = new_sequence_from_bam(alig);

	//Get cycles and positions
	cycles = alig->core.l_qseq;
	init_pos = alig->core.pos + 1;
	end_pos = alig->core.pos + cycles;

	//Duplicates check
	if(l_ult_seq)
	{
		if(pos_ult_seq == init_pos && l_ult_seq == cycles && strcmp(bam_seq, ult_seq) == 0)
		{
			//printf("\nDUPLICATE POS: %d CYCLES: %d\n\tSEQ:  %s\n\tLAST: %s", init_pos, cycles, bam_seq, ult_seq);
			duplicated++;
			_mm_free(bam_seq);
			return NO_ERROR;
		}
	}

	//Get quals
	quals = new_quality_from_bam(alig, 0);
	//bam_seq = aux_alig->sequence;
	//quals = aux_alig->quality;

	//Allocations
	ref_seq = (uint8_t *)_mm_malloc((cycles + 1) * sizeof(uint8_t), MEM_ALIG_SIZE);
	comp_res = (uint8_t *)_mm_malloc(cycles * sizeof(uint8_t), MEM_ALIG_SIZE);
	dinucs = (uint8_t *)_mm_malloc(cycles * sizeof(uint8_t), MEM_ALIG_SIZE);

	//Obtain reference for this 100 nucleotides
	flag = (uint32_t) alig->core.flag;

	genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, alig->core.tid, &init_pos, &end_pos, ref);

	//Iterates nucleotides in this read
	for(i = 0; i < cycles; i++)
	{
		/*#ifdef __SSE2__ //SSE Block
		if( (i + 16) < alig->core.l_qseq)
		{
			//Use SSE
			//_mm_prefetch(&ref_seq[i + 16], _MM_HINT_T0);
			//_mm_prefetch(&bam_seq[i + 16], _MM_HINT_T0);

			//Pack sequences
			v_ref = _mm_load_si128(&ref_seq[i]);
			v_seq = _mm_load_si128(&bam_seq[i]);

			//Compare sequences
			v_comp = _mm_cmpeq_epi8(v_ref, v_seq);

			//Store comparation values
			_mm_store_si128(&comp_res[i], v_comp);

			i += 15;
		}
		else
		#endif //SSE Block*/
		{
			if(ref_seq[i] != bam_seq[i])
			{
				comp_res[i] = 0x00;	//Diff
			}
			else
			{
				comp_res[i] = 0xFF;	//Equals
			}
		}
	}

	//Dinucs
	for(i = 0; i < cycles; i++)
	{
		if(i > 0)
		{
			recal_get_dinuc(bam_seq[i-1], bam_seq[i], &dinucs[i]);
		}
		else
		{
			dinucs[i] = d_X;
		}
	}

	//Add data
	recal_add_base_v(output_data, bam_seq, quals, 0, cycles, dinucs, comp_res);


	//Set last sequence for duplicates
	strcpy(ult_seq, bam_seq);
	l_ult_seq = cycles;
	pos_ult_seq = init_pos;

	_mm_free(ref_seq);
	_mm_free(comp_res);
	_mm_free(dinucs);
	_mm_free(bam_seq);
	free(quals);
	//alignment_free(aux_alig);


	/*FILE *fp;
	uint8_t *quals;
	uint8_t *ref_seq;
	char *bam_seq, *bam_seq_aux;
	uint8_t aux_comp[16];
	uint8_t *comp_res, *dinucs;
	uint8_t i, j;
	uint8_t match;
	#ifdef __SSE2__
	__m128i v_ref, v_seq, v_comp;
	#endif
	DINUCLEOTIDE dinuc;
	unsigned long int init_pos, end_pos;
	uint32_t flag;

	//Convert bam1_t to alignment_t
	//aux_alig = alignment_new_by_bam(alig, 1);

	//Unmapped readings
	if(alig->core.tid < 0)
	{
		unmapped++;
		return NO_ERROR;
	}

	//Bam seq fields
	bam_seq = convert_to_sequence_string(bam1_seq(alig), alig->core.l_qseq);

	//Qual fields
	quals = (uint8_t *)calloc(alig->core.l_qseq + 1, sizeof(uint8_t));
	convert_to_quality_string_length(quals, bam1_qual(alig), alig->core.l_qseq, 1);

	//Sequence fields
	ref_seq = (uint8_t *)malloc(alig->core.l_qseq * sizeof(uint8_t));
	flag = (uint32_t) alig->core.flag;
	init_pos = alig->core.pos + 1;
	end_pos = init_pos + (alig->core.l_qseq - 1);

	//Obtain reference for this 100 nucleotides
	genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, alig->core.tid, &init_pos, &end_pos, ref);
ref_seq
	//DEBUG
	#ifdef D_SAMPLES_OUTPUT
	if(rand() % D_PROB_SAMPLE == 1)
	{
		fp = fopen(out_path, "a+");
		fprintf(fp, "--------------------------\n");
		fprintf(fp, "Chrom:%d\nPos ref:%lu\nPos bam:%d\nL_qseq:%d\nSEQ:%s\nBAM:%s\n", alig->core.tid, init_pos, alig->core.pos, alig->core.l_qseq, ref_seq, bam_seq);
		fprintf(fp, "Mis:");
		for(i = 0; i < alig->core.l_qseq; i++)
		{
			fprintf(fp, "%d",  ref_seq[i] != bam_seq[i]);
		}
		fprintf(fp,"\n");
		fclose(fp);
	}
	#endif

	//Iterates nucleotides in this read
	dinuc = 0;
	comp_res = malloc(alig->core.l_qseq * sizeof(uint8_t));
	dinucs = malloc(alig->core.l_qseq * sizeof(uint8_t));
	for(i = 0; i < alig->core.l_qseq; i++)
	{
		#ifdef __SSE2adfas__ //SSE Block
		if( (i + 16) < alig->core.l_qseq)
		{
			//Use SSE
			//_mm_prefetch(&ref_seq[i + 16], _MM_HINT_T0);
			//_mm_prefetch(&bam_seq[i + 16], _MM_HINT_T0);

			//Pack sequences
			v_ref = _mm_load_si128(&ref_seq[i]);
			v_seq = _mm_load_si128(&bam_seq[i]);

			//Compare sequences
			v_comp = _mm_cmpeq_epi8(v_ref, v_seq);

			//Store comparation values
			_mm_store_si128(&comp_res[i], v_comp);

			for(j = i; j < i + 16; j++)
			{
				#ifdef NOT_COUNT_NUCLEOTIDE_N
				if(bam_seq[j] != 'N')
				#endif
				{
					//Only take dinuc if we arent in first cycle
					if(j > 0)
					{
						dinuc = recal_get_dinuc(bam_seq[j-1], bam_seq[j]);
					}
					else
					{
						dinuc = recal_get_dinuc('_', '_');
					}

					//if this happen, then error will be thrown
					if(dinuc == -1)
					{
						dinuc = d_X;
						//printf("\nDinuc:%c%c", bam_seq[i-1],bam_seq[i]);
					}

					//Add data
					recal_add_base(output_data, quals[j], j, dinuc, v_comp[j-i]);
				}
			}
			i += 15;
		}
		else
		#endif //SSE Block
		{
			#ifdef NOT_COUNT_NUCLEOTIDE_N
			if(bam_seq[i] != 'N')
			#endif
			{
				//Normal
				//Only take dinuc if we arent in first cycle
				if(i > 0)
				{
					dinuc = recal_get_dinuc(bam_seq[i-1], bam_seq[i]);
				}
				else
				{
					dinuc = recal_get_dinuc('_', '_');
				}

				//if this happen, then error will be thrown
				if(dinuc == -1)
				{
					dinuc = d_X;
					//printf("\nDinuc:%c%c", bam_seq[i-1],bam_seq[i]);
				}

				//Add data
				recal_add_base(output_data, quals[i], i, dinuc, ref_seq[i] != bam_seq[i]);
			}
			if(ref_seq[i] != bam_seq[i])
			{
				comp_res[i] = 0x00;	//Diff
			}
			else
			{
				comp_res[i] = 0xFF;	//Equals
			}
		}
	}

	//Dinucs
	for(i = 0; i < alig->core.l_qseq; i++)
	{
		if(i > 0)
		{
			recal_get_dinuc(bam_seq[i-1], bam_seq[i], &dinucs[i]);
			if(dinucs[i] == -1)
			{
				dinucs[i] = d_X;
			}
		}
		else
		{
			dinucs[i] = d_X;
		}
	}

	//Add data
	recal_add_base_v(output_data, bam_seq, quals, 0, alig->core.l_qseq - 1, dinucs, comp_res);

	free(comp_res);
	free(dinucs);
	free(quals);
	free(ref_seq);
	free(bam_seq);*/

	//alignment_free(aux_alig);

	return NO_ERROR;
}
