#include "data_collector.h"

/**
 * Get recalibration data from BAM path.
 */
void recal_get_data_from_file(char *bam_path, char *ref_name, char *ref_path, recal_info_t *out_info)
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
		time_init_slot(D_SLOT_GET_DATA_BAM, clock(), time_global_stats);
	#endif
	recal_get_data_from_bam(bam_f, ref, out_info);
	#ifdef D_TIME_DEBUG
		time_set_slot(D_SLOT_GET_DATA_BAM, clock(), time_global_stats);
	#endif

	//Memory free
	printf("\nClosing BAM file...\n");
	bam_fclose(bam_f);
	genome_free(ref);
	printf("BAM closed.\n");
}

/**
 * Get recalibration data from BAM file.
 */
void recal_get_data_from_bam(bam_file_t *bam, genome_t* ref, recal_info_t* output_data)
{
	bam_batch_t* batch;
	int count = 0;

	#ifdef USE_BATCH_POOL
	bam_pool_t* pool;
	#endif

	int zero = 0;

	unmapped = 0;

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
		time_init_slot(D_SLOT_READ_BATCH, clock(), time_global_stats);
	#endif
	bam_fread_max_size(batch, MAX_BATCH_SIZE, 1, bam);
	#ifdef D_TIME_DEBUG
		time_set_slot(D_SLOT_READ_BATCH, clock(), time_global_stats);
	#endif

	printf("\nNum alignments in batchs: %d\n----------------\n", batch->num_alignments);

	#ifdef USE_BATCH_POOL
	printf("Using batch pool!\n");
	pool = pool_new(3, MAX_BATCH_SIZE);
	#endif

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
			time_init_slot(D_SLOT_GET_DATA_BATCH, clock(), time_global_stats);
		#endif
		recal_get_data_from_bam_batch(batch, ref, output_data);
		#ifdef D_TIME_DEBUG
			time_set_slot(D_SLOT_GET_DATA_BATCH, clock(), time_global_stats);
		#endif

		//Update read counter
		count += batch->num_alignments;

		//Show total progress
		printf("\n Total alignments readed: %d", count);

		//Free memory and take a new batch
		#ifndef USE_BATCH_POOL
		bam_batch_free(batch, 1);
		#endif
		batch = bam_batch_new(MAX_BATCH_SIZE, MULTIPLE_CHROM_BATCH);

		//Read batch
		#ifdef D_TIME_DEBUG
			time_init_slot(D_SLOT_READ_BATCH, clock(), time_global_stats);
		#endif
		bam_fread_max_size(batch, MAX_BATCH_SIZE, 1, bam);
		#ifdef D_TIME_DEBUG
			time_set_slot(D_SLOT_READ_BATCH, clock(), time_global_stats);
		#endif
	}

	printf("\n----------------\n%d alignments readed.", count);
	printf("\n%d alignments processed.", count - unmapped);
	printf("\n%d alignments unmapped.", unmapped);

	#ifdef USE_BATCH_POOL
		pool_free(pool);
	#endif

	//Free batch
	bam_batch_free(batch, 1);
}

/**
 * Get recalibration data from BAM batch of alignments.
 */
void recal_get_data_from_bam_batch(bam_batch_t* batch, genome_t* ref, recal_info_t* output_data)
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
			time_init_slot(D_SLOT_GET_DATA_ALIG, clock(), time_global_stats);
		#endif
		recal_get_data_from_bam_alignment(batch->alignments_p[i], ref, output_data);
		#ifdef D_TIME_DEBUG
			time_set_slot(D_SLOT_GET_DATA_ALIG, clock(), time_global_stats);
		#endif
		}
	}
}

/**
 * Get recalibration data from alignment.
 */
void recal_get_data_from_bam_alignment(bam1_t* alig, genome_t* ref, recal_info_t* output_data)
{
	//alignment_t* aux_alig;
	FILE *fp;
	char *quals;
	char *ref_seq;
	char *bam_seq, *bam_seq_aux;
	char aux_comp[16];
	char *comp_res, *dinucs;
	int i, miss, j;
	#ifdef __SSE2__
	__m128i v_ref, v_seq, v_comp;
	#endif
	enum DINUC dinuc;
	unsigned long int init_pos, end_pos;
	uint32_t flag;

	//Convert bam1_t to alignment_t
	//aux_alig = alignment_new_by_bam(alig, 1);

	//Unmapped readings
	if(alig->core.tid < 0)
	{
		unmapped++;
		return;
	}

	//Bam seq fields
	bam_seq = convert_to_sequence_string(bam1_seq(alig), alig->core.l_qseq);

	//Qual fields
	quals = (char*) calloc(alig->core.l_qseq + 1, sizeof(char));
	convert_to_quality_string_length(quals, bam1_qual(alig), alig->core.l_qseq, 1);

	//Sequence fields
	ref_seq = malloc(alig->core.l_qseq * sizeof(char));
	flag = (uint32_t) alig->core.flag;
	init_pos = alig->core.pos + 1;
	end_pos = init_pos + (alig->core.l_qseq - 1);

	//Obtain reference for this 100 nucleotides
	genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, alig->core.tid, &init_pos, &end_pos, ref);

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
	comp_res = malloc(alig->core.l_qseq * sizeof(char));
	dinucs = malloc(alig->core.l_qseq * sizeof(int));
	for(i = 0; i < alig->core.l_qseq; i++)
	{
		#ifdef __SSE2__ /*SSE Block*/
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

			/*for(j = i; j < i + 16; j++)
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
			}*/
			i += 15;
		}
		else
		#endif /*SSE Block*/
		{
			/*#ifdef NOT_COUNT_NUCLEOTIDE_N
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
			}*/
			if(ref_seq[i] != bam_seq[i])
			{
				comp_res[i] = 0;	/*Diff*/
			}
			else
			{
				comp_res[i] = -1;	/*Equals*/
			}
		}
	}

	//Dinucs
	for(i = 0; i < alig->core.l_qseq; i++)
	{
		if(i > 0)
		{
			dinucs[i] = recal_get_dinuc(bam_seq[i-1], bam_seq[i]);
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
	free(bam_seq);

	//alignment_free(aux_alig);
}
