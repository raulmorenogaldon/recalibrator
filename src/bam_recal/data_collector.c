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
	#ifdef NOT_PRIMARY_ALIGNMENT
		notprimary = 0;
	#endif
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
	printf("\n%d alignments processed.", count - unmapped - mapzero - duplicated - notprimary);
	printf("\n%d alignments duplicated.", duplicated);
	#ifdef NOT_PRIMARY_ALIGNMENT
		printf("\n%d not primary alignments.", notprimary);
	#endif
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
		#ifdef D_TIME_DEBUG
			time_init_slot(D_SLOT_GET_DATA_ALIG, clock(), TIME_GLOBAL_STATS);
		#endif
		//Recollection
		recal_get_data_from_bam_alignment(batch->alignments_p[i], ref, output_data);
		#ifdef D_TIME_DEBUG
			time_set_slot(D_SLOT_GET_DATA_ALIG, clock(), TIME_GLOBAL_STATS);
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
	char *ref_seq;
	uint8_t aux_comp[16];
	size_t init_pos, end_pos;
	uint32_t cycles;
	uint8_t *comp_res;
	uint8_t *dinucs;
	uint32_t flag;

	//Cigar
	char *cigar;

	#ifdef __SSE2__
	__m128i v_ref, v_seq, v_comp;
	#endif

	unsigned int i, j;

	//FILTERS
	{
		//Map quality is zero
		#ifdef NOT_MAPPING_QUAL_ZERO
		if(alig->core.qual == 0)
		{
			mapzero++;
			return NO_ERROR;
		}
		#endif

		//Not primary
		#ifdef NOT_PRIMARY_ALIGNMENT
		if(alig->core.flag & 256)	//Not primary alignment flag
		{
			notprimary++;
			return NO_ERROR;
		}
		#endif

		//Unmapped readings
		if(alig->core.tid < 0)
		{
			unmapped++;
			return NO_ERROR;
		}
	}

	//aux_alig = alignment_new_by_bam(alig, 1);

	//Get sequence
	bam_seq = new_sequence_from_bam(alig);

	//Decompose cigar for indels
	//TODO
	//cigar = alig->core
	//decompose_cigar(, uint8_t cigar_l, char *n_elem, char *type, uint8_t *types_l, uint8_t max_types_length);

	//Get cycles and positions
	cycles = alig->core.l_qseq;
	init_pos = alig->core.pos + 1;
	end_pos = alig->core.pos + cycles;

	//Duplicates check
	/*if(l_ult_seq)
	{
		if(pos_ult_seq == init_pos && l_ult_seq == cycles && strcmp(bam_seq, ult_seq) == 0)
		{
			//printf("\nDUPLICATE POS: %d CYCLES: %d\n\tSEQ:  %s\n\tLAST: %s", init_pos, cycles, bam_seq, ult_seq);
			duplicated++;
			_mm_free(bam_seq);
			return NO_ERROR;
		}
	}*/

	//Get quals
	quals = new_quality_from_bam(alig, 0);
	//bam_seq = aux_alig->sequence;
	//quals = aux_alig->quality;

	//Allocations
	ref_seq = (char *)_mm_malloc((cycles + 1) * sizeof(char), MEM_ALIG_SIZE);
	comp_res = (uint8_t *)_mm_malloc(cycles * sizeof(uint8_t), MEM_ALIG_SIZE);
	dinucs = (uint8_t *)_mm_malloc(cycles * sizeof(uint8_t), MEM_ALIG_SIZE);

	//Obtain reference for this 100 nucleotides
	flag = (uint32_t) alig->core.flag;

	genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, (unsigned int)alig->core.tid, &init_pos, &end_pos, ref);

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

	return NO_ERROR;
}
