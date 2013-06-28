#include "bam_recal.h"

/**
 * Recalibrate BAM file from path and store in file.
 */
void recal_recalibrate_bam_file(char *orig_bam_path, recal_info_t *bam_info, char *recal_bam_path)
{
	bam_file_t *orig_bam_f, *recal_bam_f;
	bam_header_t *recal_bam_header;

	//Open bam
	printf("Opening BAM from \"%s\" to being recalibrated ...\n", orig_bam_path);
	orig_bam_f = bam_fopen(orig_bam_path);
	printf("BAM opened!...\n");

	//Create new bam
	printf("Creating new bam file in \"%s\"...\n", recal_bam_path);
	recal_bam_header = create_empty_bam_header(orig_bam_f->bam_header_p->n_targets);
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
}

/**
 * Recalibrate BAM file and store in file.
 */
void recal_recalibrate_bam(bam_file_t *orig_bam_f, recal_info_t *bam_info, bam_file_t *recal_bam_f)
{
	bam_batch_t* batch;
	int count = 0;

	#ifdef USE_BATCH_POOL
	bam_pool_t* pool;
	#endif

	//Allocate memory for batchs
	batch = bam_batch_new(MAX_BATCH_SIZE, MULTIPLE_CHROM_BATCH);

	//Read batch
	bam_fread_max_size(batch, MAX_BATCH_SIZE, 1, orig_bam_f);

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
		bam_fread_max_size(batch, MAX_BATCH_SIZE, 1, orig_bam_f);

	}
	#ifdef D_RECAL_INTER_RESULTS
		printf("\n=======================\n");
	#endif

	#ifdef USE_BATCH_POOL
		pool_free(pool);
	#endif

	printf("\n---------------------\n", count);
}

/**
 * Recalibrate BAM batch of alignments and store in file.
 */
void recal_recalibrate_batch(bam_batch_t* batch, recal_info_t *bam_info, bam_file_t *recal_bam_f)
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
}

/**
 * Recalibrate alignment and store in file.
 */
void recal_recalibrate_alignment(bam1_t* alig, recal_info_t *bam_info, bam_file_t *recal_bam_f)
{
	int qual_index;
	int matrix_index;
	int i;
	enum DINUC dinuc;
	double delta_r, delta_rc, delta_rd;
	char *quals;
	char *bam_seq;

	alignment_t* aux_alig;
	bam1_t *aux_alig1;

	//Convert bam1_t to alignment_t
	aux_alig = alignment_new_by_bam(alig, 1);

	//Bam seq fields
	bam_seq = convert_to_sequence_string(bam1_seq(alig), alig->core.l_qseq);

	//Qual fields
	quals = (char*) calloc(alig->core.l_qseq + 1, sizeof(char));
	convert_to_quality_string_length(quals, bam1_qual(alig), alig->core.l_qseq, 1);

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
				dinuc = recal_get_dinuc(bam_seq[i-1], bam_seq[i]);
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
			quals[i] = quals[i] + bam_info->total_delta + delta_r + delta_rc + delta_rd;

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

	//Set qualities in alignment
	free(aux_alig->quality);
	aux_alig->quality = quals;
	aux_alig1 = convert_to_bam(aux_alig, 1);

	//Write in bam
	bam_fwrite(aux_alig1, recal_bam_f);

	//Memory free
	alignment_free(aux_alig);
	bam_destroy1(aux_alig1);
	free(bam_seq);
}
