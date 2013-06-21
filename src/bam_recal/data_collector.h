#ifndef DATA_COLLECTOR_H_
#define DATA_COLLECTOR_H_

#include "genome.h"
#include "alignment.h"
#include "bam_file.h"
#include "config.h"

#ifdef D_TIME_DEBUG
	#include "timestats.h"
#endif

/**
 * Get recalibration data from BAM path.
 */
void recal_get_data_from_file(char *bam_path, char *ref_name, char *ref_path, recal_info_t *out_info);

/**
 * Get recalibration data from BAM file.
 */
void recal_get_data_from_bam(bam_file_t *bam, genome_t* ref, recal_info_t* output_data);

/**
 * Get recalibration data from BAM batch of alignments.
 */
void recal_get_data_from_bam_batch(bam_batch_t* batch, genome_t* ref, recal_info_t* output_data);

/**
 * Get recalibration data from alignment.
 */
void recal_get_data_from_bam_alignment(bam1_t* alig, genome_t* ref, recal_info_t* output_data);

#endif /* DATA_COLLECTOR_H_ */
