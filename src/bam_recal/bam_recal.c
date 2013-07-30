#include "bam_recal.h"

/**
 * Recalibrate BAM file from path and store in file.
 */
ERROR_CODE
recal_recalibrate_bam_file(const char *orig_bam_path, const recal_info_t *bam_info, const char *recal_bam_path)
{
	bam_file_t *orig_bam_f, *recal_bam_f;
	bam_header_t *recal_bam_header;

	//Open bam
	printf("Opening BAM from \"%s\" to being recalibrated ...\n", orig_bam_path);
	orig_bam_f = bam_fopen(orig_bam_path);
	printf("BAM opened!...\n");

	//Allocate
	recal_bam_header = (bam_header_t *) calloc(1, sizeof(bam_header_t));

	//Create new bam
	printf("Creating new bam file in \"%s\"...\n", recal_bam_path);
	init_empty_bam_header(orig_bam_f->bam_header_p->n_targets, recal_bam_header);
	recal_bam_f = bam_fopen_mode(recal_bam_path, recal_bam_header, "w");
	bam_fwrite_header(recal_bam_header, recal_bam_f);
	printf("New BAM initialized!...\n");

	//Recalibrate bams
	#ifdef D_TIME_DEBUG
		time_init_slot(D_SLOT_RECALIBRATE, clock(), TIME_GLOBAL_STATS);
	#endif
	recal_recalibrate_bam(orig_bam_f, bam_info, recal_bam_f);
	#ifdef D_TIME_DEBUG
		time_set_slot(D_SLOT_RECALIBRATE, clock(), TIME_GLOBAL_STATS);
	#endif

	//Memory free
	printf("Closing \"%s\" BAM file...\n", orig_bam_path);
	bam_fclose(orig_bam_f);
	printf("Closing \"%s\" BAM file...\n", recal_bam_path);
	bam_fclose(recal_bam_f);
	printf("BAMs closed.\n");
	printf("Recalibration DONE.\n");

	return NO_ERROR;
}

/**
 * Recalibrate BAM file and store in file.
 */
ERROR_CODE
recal_recalibrate_bam(const bam_file_t *orig_bam_f, const recal_info_t *bam_info, bam_file_t *recal_bam_f)
{
	bam_batch_t* batch;
	int count = 0;

	#ifdef USE_BATCH_POOL
	bam_pool_t* pool;
	#endif

	//Allocate memory for batchs
	batch = bam_batch_new(MAX_BATCH_SIZE, MULTIPLE_CHROM_BATCH);

	//Read batch
	bam_fread_max_size(batch, MAX_BATCH_SIZE, 0, orig_bam_f);

	#ifdef USE_BATCH_POOL
	printf("Using batch pool!\n");
	pool = pool_new(3, MAX_BATCH_SIZE);
	#endif

	#ifdef D_RECAL_INTER_RESULTS
		printf("=======================\nQual intermediate results\n");
	#endif

	while(batch->num_alignments != 0
		#ifdef D_MAX_READS_W
			&& count < D_MAX_READS_W
		#endif
		)
	{
		#ifdef USE_BATCH_POOL
		if(pool_has_free_batchs(pool))
		{
			pool_add_batch(batch, pool);
			bam_batch_free(batch, 1);

			//Get batch in pool
			pool_next_batch(pool);
			batch = pool_get_batch(pool);
		}
		#endif

		//Recalibrate batch
		recal_recalibrate_batch(batch, bam_info, recal_bam_f);
		//Update read counter
		count += batch->num_alignments;

		//Show total progress
		printf("\n\tTotal alignments recalibrated: %d\n", count);

		//Free memory and take a new batch
		#ifndef USE_BATCH_POOL
		bam_batch_free(batch, 1);
		#endif
		batch = bam_batch_new(MAX_BATCH_SIZE, MULTIPLE_CHROM_BATCH);

		//Read batch
		bam_fread_max_size(batch, MAX_BATCH_SIZE, 0, orig_bam_f);

	}
	#ifdef D_RECAL_INTER_RESULTS
		printf("\n=======================\n");
	#endif

	#ifdef USE_BATCH_POOL
		pool_free(pool);
	#endif

	printf("\n---------------------\n", count);

	return NO_ERROR;
}

/**
 * Recalibrate BAM batch of alignments and store in file.
 */
ERROR_CODE
recal_recalibrate_batch(const bam_batch_t* batch, const recal_info_t *bam_info, bam_file_t *recal_bam_f)
{
	int i;

	//Process all alignments of the batchs
	for(i = 0; i < batch->num_alignments; i++)
	{
		#ifdef D_TIME_DEBUG
		time_init_slot(D_SLOT_RECAL_ALIG, clock(), TIME_GLOBAL_STATS);
		#endif
		//Process every alignment
		recal_recalibrate_alignment(batch->alignments_p[i], bam_info, recal_bam_f);
		#ifdef D_TIME_DEBUG
		time_set_slot(D_SLOT_RECAL_ALIG, clock(), TIME_GLOBAL_STATS);
		#endif
	}

	return NO_ERROR;
}

/**
 * Recalibrate alignment and store in file.
 */
ERROR_CODE
recal_recalibrate_alignment(const bam1_t* alig, const recal_info_t *bam_info, bam_file_t *recal_bam_f)
{
	unsigned int qual_index;
	unsigned int matrix_index;
	unsigned int i;
	uint8_t dinuc;
	double delta_r, delta_rc, delta_rd;
	char *quals, *res_quals;
	char *bam_seq;
	double estimated_Q;

	alignment_t* aux_alig;
	bam1_t *aux_alig1;

	//Bam seq fields
	bam_seq = new_sequence_from_bam(alig);

	//Qual fields
	quals = new_quality_from_bam(alig, 0);
	res_quals = (char *)malloc(alig->core.l_qseq * sizeof(char));
	//quals = (char*) calloc(alig->core.l_qseq + 1, sizeof(char));
	//convert_to_quality_string_length(quals, bam1_qual(alig), alig->core.l_qseq, 1);

	//Iterates nucleotides in this read
	dinuc = 0;
	for(i = 0; i < alig->core.l_qseq; i++)
	{
		//Compare only if the nucleotide is not "N"
		#ifdef NOT_COUNT_NUCLEOTIDE_N
		if(bam_seq[i] != 'N')
		#endif
		{
			//Recalibrate quality
			qual_index = quals[i] - bam_info->min_qual;
			delta_r = bam_info->qual_delta[qual_index];

			matrix_index = qual_index * bam_info->num_cycles + i;
			delta_rc = bam_info->qual_cycle_delta[matrix_index];

			//dont take prev dinuc in first cycle (delta = 0)
			if(i > 0)
			{
				recal_get_dinuc(bam_seq[i-1], bam_seq[i], &dinuc);
				matrix_index = qual_index * bam_info->num_dinuc + i;
				delta_rd = bam_info->qual_dinuc_delta[matrix_index];
			}
			else
			{
				delta_rd = 0.0;
			}

			#ifdef D_RECAL_INTER_RESULTS
				int aux_qual = quals[i];
			#endif

			//Recalibration formula
				double global_delta = bam_info->total_delta;
				double calidad = (double)quals[i];
				if(calidad > MIN_QUALITY_TO_STAT)
				{
					//double res = (double)quals[i] + bam_info->total_delta + delta_r + delta_rc + delta_rd;
					recal_get_estimated_Q(bam_info->qual_bases, bam_info->num_quals, 0, &estimated_Q);
					double res = estimated_Q + bam_info->total_delta + delta_r + delta_rc + delta_rd;
					res_quals[i] = (char)res;
				}
				else
				{
					res_quals[i] = calidad;
				}
			//quals[i] = (char)((double)quals[i] + bam_info->total_delta + delta_r + delta_rc + delta_rd);

			#ifdef D_RECAL_INTER_RESULTS
			if(rand() % D_RECAL_INTER_PROB == 1)
			{
				printf("Q: %d ",aux_qual);
				printf("---- %d\n",quals[i]);
				printf("Cy: %d, Qr + DQg:%f + DQr:%f + DQc:%f + DQd:%f\n", i, bam_info->total_delta, delta_r, delta_rc, delta_rd);
			}
			#endif
		}
	}

	//Convert bam1_t to alignment_t
	aux_alig = alignment_new_by_bam(alig, 0);

	//Set qualities in alignment
	free(aux_alig->quality);
	aux_alig->quality = res_quals;
	//Fix reads in alignment (sequence conversion to string is bug)
	free(aux_alig->sequence);
	aux_alig->sequence = bam_seq;

	aux_alig1 = convert_to_bam(aux_alig, 0);
	bam_seq = new_sequence_from_bam(aux_alig1);
	quals = new_quality_from_bam(aux_alig1, 0);
	//Write in bam
	bam_fwrite(aux_alig1, recal_bam_f);

	//Memory free
	free(quals);
	alignment_free(aux_alig);
	bam_destroy1(aux_alig1);

	return NO_ERROR;
}
