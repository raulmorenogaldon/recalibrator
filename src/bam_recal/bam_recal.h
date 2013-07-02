#ifndef BAM_RECAL_H_
#define BAM_RECAL_H_

#include <bam_file.h>
#include <recal_common.h>
#include <recal_config.h>
#include <aux_library.h>
#include <timestats.h>

#include "recal_structs.h"

/***********************************************
 * BAM RECALIBRATION PHASE 2 - RECALIBRATION
 **********************************************/

/**
 * Recalibrate BAM file from path and store in file.
 */
extern void recal_recalibrate_bam_file(char *orig_bam_path, recal_info_t *bam_info, char *recal_bam_path);

/**
 * Recalibrate BAM file and store in file.
 */
extern void recal_recalibrate_bam(bam_file_t *orig_bam_f, recal_info_t *bam_info, bam_file_t *recal_bam_f);

/**
 * Recalibrate BAM batch of alignments and store in file.
 */
extern void recal_recalibrate_batch(bam_batch_t* batch, recal_info_t *bam_info, bam_file_t *recal_bam_f);

/**
 * Recalibrate alignment and store in file.
 */
extern void recal_recalibrate_alignment(bam1_t* alig, recal_info_t *bam_info, bam_file_t *recal_bam_f);

#endif /* BAM_RECAL_H_ */
