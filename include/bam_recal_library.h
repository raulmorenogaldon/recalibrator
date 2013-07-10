/*
 * bam_recal_library.h
 *
 *  Created on: Jun 28, 2013
 *      Author: rmoreno
 */

#ifndef BAM_RECAL_LIBRARY_H_
#define BAM_RECAL_LIBRARY_H_

#include <bam.h>
#include <bam_file.h>
#include <alignment.h>
#include <genome.h>
#include <recal_common.h>

/**
 * Recalibration data storage
 * This struct hold all data necessary to perform recalibrations
 */
struct recal_info;
typedef struct recal_info recal_info_t;

/**
 * Dinucleotide enumeration.
 * Represents all possible dinucleotide cases.
 */
typedef enum DINUC
{
	dAA = 0,
	dAG = 1,
	dAC = 2,
	dAT = 3,
	dGA = 4,
	dGG = 5,
	dGC = 6,
	dGT = 7,
	dCA = 8,
	dCG = 9,
	dCC = 10,
	dCT = 11,
	dTA = 12,
	dTG = 13,
	dTC = 14,
	dTT = 15,
	d_X = 16	/**< Not a dinucleotide. For example "NT" or "NN".*/
} DINUCLEOTIDE;


/***********************************************
 * DATA MANAGEMENT
 **********************************************/

/**
 * Allocate new recalibration data.
 */
EXTERNC ERROR_CODE recal_new_info(const int cycles, recal_info_t **out_info);

/**
 * Free recalibration data.
 */
EXTERNC ERROR_CODE recal_destroy_info(recal_info_t **data);

/**
 * Add recalibration data from one base.
 */
EXTERNC ERROR_CODE recal_add_base(recal_info_t *data, const int qual, const int cycle, const int dinuc, const int miss);

/**
 * Add recalibration data from vector of bases
 */
EXTERNC ERROR_CODE recal_add_base_v(recal_info_t *data, const char *seq, const char *quals, const int init_cycle, const int end_cycle, const char *dinuc, const char *misses);

/**
 * Compute deltas from bases and misses.
 */
EXTERNC ERROR_CODE recal_calc_deltas(recal_info_t* data);


/***********************************************
 * DINUC ENUM OPERATIONS
 **********************************************/

/**
 * Return dinucleotide enumeration from two bases.
 */
EXTERNC ERROR_CODE recal_get_dinuc(const char A, const char B, DINUCLEOTIDE *out_dinuc);


/***********************************************
 * BAM RECALIBRATION PHASE 1 - DATA COLLECT
 **********************************************/

/**
 * Get recalibration data from BAM path.
 */
EXTERNC ERROR_CODE recal_get_data_from_file(char *bam_path, char *ref_name, char *ref_path, recal_info_t *out_info);

/**
 * Get recalibration data from BAM file.
 */
EXTERNC ERROR_CODE recal_get_data_from_bam(bam_file_t *bam, genome_t* ref, recal_info_t* output_data);

/**
 * Get recalibration data from BAM batch of alignments.
 */
EXTERNC ERROR_CODE recal_get_data_from_bam_batch(bam_batch_t* batch, genome_t* ref, recal_info_t* output_data);

/**
 * Get recalibration data from alignment.
 */
EXTERNC ERROR_CODE recal_get_data_from_bam_alignment(bam1_t* alig, genome_t* ref, recal_info_t* output_data);


/***********************************************
 * BAM RECALIBRATION PHASE 2 - RECALIBRATION
 **********************************************/

/**
 * Recalibrate BAM file from path and store in file.
 */
EXTERNC ERROR_CODE recal_recalibrate_bam_file(char *orig_bam_path, recal_info_t *bam_info, char *recal_bam_path);

/**
 * Recalibrate BAM file and store in file.
 */
EXTERNC ERROR_CODE recal_recalibrate_bam(bam_file_t *orig_bam_f, recal_info_t *bam_info, bam_file_t *recal_bam_f);

/**
 * Recalibrate BAM batch of alignments and store in file.
 */
EXTERNC ERROR_CODE recal_recalibrate_batch(bam_batch_t* batch, recal_info_t *bam_info, bam_file_t *recal_bam_f);

/**
 * Recalibrate alignment and store in file.
 */
EXTERNC ERROR_CODE recal_recalibrate_alignment(bam1_t* alig, recal_info_t *bam_info, bam_file_t *recal_bam_f);


/***********************************************
 * FILE OPERATIONS
 **********************************************/

/**
 * Print to file data from recalibration.
 */
EXTERNC ERROR_CODE recal_fprint_info(const recal_info_t *data, const char *path);

/**
 * Save to file recalibration data.
 */
EXTERNC ERROR_CODE recal_save_recal_info(const recal_info_t *data, const char *path);

/**
 * Load from file recalibration data.
 */
EXTERNC ERROR_CODE recal_load_recal_info(const char *path, recal_info_t *data);

#endif /* BAM_RECAL_LIBRARY_H_ */
