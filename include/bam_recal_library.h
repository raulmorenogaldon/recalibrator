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
 * \brief Recalibration data storage.
 *
 * This struct hold all data necessary to perform recalibrations.
 */
struct recal_info;
typedef struct recal_info recal_info_t;

/**
 * \brief Dinucleotide enumeration.
 *
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
	d_X = 16	/* Not a dinucleotide. For example "NT" or "NN".*/
} DINUCLEOTIDE;


/***********************************************
 * DATA MANAGEMENT
 **********************************************/

/**
 * \brief Initialize empty recalibration data struct.
 *
 * \param cycles Number of maximum cycles to stat.
 * \param out_info Previously allocated info struct to initialize.
 */
EXTERNC ERROR_CODE recal_init_info(const uint32_t cycles, recal_info_t **out_data);

/**
 * \brief Free all resources of recalibration.
 *
 * Free all data struct attributes including itself at the end.
 *
 * \param data Data struct to free
 */
EXTERNC ERROR_CODE recal_destroy_info(recal_info_t **data);

/**
 * \brief Add recalibration data from one base.
 *
 * \param data Data struct to add stats.
 * \param qual Quality to add.
 * \param cycle Cycle to add.
 * \param dinuc Dinucleotide to add.
 * \param miss Indicate if match(!=0) or not (0)
 */
EXTERNC ERROR_CODE recal_add_base(recal_info_t *data, const uint8_t qual, const uint16_t cycle, const uint8_t dinuc, const double match) __ATTR_HOT;

/**
 * \brief Add recalibration data from vector of bases.
 *
 * \param data Vector of data structs to add stats.
 * \param qual Vector of qualities to add.
 * \param cycle Init cycle to add.
 * \param cycle End cycles to add.
 * \param dinuc Vector of dinucleotides to add.
 * \param miss Vector of match(!=0) or not (0)
 */
EXTERNC ERROR_CODE recal_add_base_v(recal_info_t *data, const uint8_t *seq, const uint8_t *quals, const uint16_t init_cycle, const uint32_t num_cycles, const uint8_t *dinuc, const double *matches) __ATTR_HOT;

/**
 * \brief Compute deltas from bases and misses.
 *
 * \param data Data to compute deltas.
 */
EXTERNC ERROR_CODE recal_calc_deltas(recal_info_t* data);


/***********************************************
 * DINUC ENUM OPERATIONS
 **********************************************/

/**
 * \brief Return dinucleotide enumeration from two bases.
 *
 * \param A First base.
 * \param B Second base.
 * \param out_dinuc Pointer to output dinucleotide.
 */
EXTERNC ERROR_CODE recal_get_dinuc(const char A, const char B, uint8_t *out_dinuc) __ATTR_HOT;


/***********************************************
 * BAM RECALIBRATION PHASE 1 - DATA COLLECT
 **********************************************/

/**
 * \brief Get recalibration data from BAM path.
 *
 * \param bam_path Path to BAM.
 * \param ref_name String with the name of the reference genome.
 * \param ref_path Path to reference.
 * \param out_info Data struct to fill.
 */
EXTERNC ERROR_CODE recal_get_data_from_file(const char *bam_path, const char *ref_name, const char *ref_path, recal_info_t *out_info);

/**
 * \brief Get recalibration data from BAM file.
 *
 * \param bam BAM file struct to process.
 * \param ref Reference genome struct.
 * \param out_info Data struct to fill.
 */
EXTERNC ERROR_CODE recal_get_data_from_bam(const bam_file_t *bam, const genome_t* ref, recal_info_t* output_data);

/**
 * \brief Get recalibration data from BAM batch of alignments.
 *
 * \param batch BAM batch struct to process.
 * \param ref Reference genome struct.
 * \param out_info Data struct to fill.
 */
EXTERNC ERROR_CODE recal_get_data_from_bam_batch(const bam_batch_t* batch, const genome_t* ref, recal_info_t* output_data);

/**
 * \brief Get recalibration data from alignment.
 *
 * \param batch BAM alignment struct to process.
 * \param ref Reference genome struct.
 * \param out_info Data struct to fill.
 */
EXTERNC ERROR_CODE recal_get_data_from_bam_alignment(const bam1_t* alig, const genome_t* ref, recal_info_t* output_data) __ATTR_HOT;


/***********************************************
 * BAM RECALIBRATION PHASE 2 - RECALIBRATION
 **********************************************/

/**
 * \brief Recalibrate BAM file from path and store in file.
 *
 * \param orig_bam_path Path to BAM which will be recalibrated.
 * \param bam_info Data struct with recalibration info.
 * \param recal_bam_path Path to output BAM.
 */
EXTERNC ERROR_CODE recal_recalibrate_bam_file(const char *orig_bam_path, const recal_info_t *bam_info, const char *recal_bam_path);

/**
 * \brief Recalibrate BAM file and store in file.
 *
 * \param orig_bam_f BAM file struct to recalibrate.
 * \param bam_info Data struct with recalibration info.
 * \param recal_bam_f Recalibrated BAM output file struct.
 */
EXTERNC ERROR_CODE recal_recalibrate_bam(const bam_file_t *orig_bam_f, const recal_info_t *bam_info, bam_file_t *recal_bam_f);

/**
 * \brief Recalibrate BAM batch of alignments and store in file.
 *
 * \param batch Batch struct to recalibrate.
 * \param bam_info Data struct with recalibration info.
 * \param recal_bam_f Recalibrated BAM output file struct.
 */
EXTERNC ERROR_CODE recal_recalibrate_batch(const bam_batch_t* batch, const recal_info_t *bam_info, bam_file_t *recal_bam_f);

/**
 * \brief Recalibrate alignment and store in file.
 *
 * \param alig BAM alignment struct to recalibrate.
 * \param bam_info Data struct with recalibration info.
 * \param recal_bam_f Recalibrated BAM output file struct.
 */
EXTERNC ERROR_CODE recal_recalibrate_alignment(const bam1_t* alig, const recal_info_t *bam_info, bam_file_t *recal_bam_f) __ATTR_HOT;


/***********************************************
 * FILE OPERATIONS
 **********************************************/

/**
 * \brief Print to file data from recalibration.
 *
 * \param data Data struct with recalibration info.
 * \param path File path to print the info.
 */
EXTERNC ERROR_CODE recal_fprint_info(const recal_info_t *data, const char *path);

/**
 * \brief Save to file recalibration data.
 *
 * \param data Data struct with recalibration info.
 * \param path File path to save data.
 */
EXTERNC ERROR_CODE recal_save_recal_info(const recal_info_t *data, const char *path);

/**
 * \brief Load from file recalibration data.
 *
 * \param path File path to load data.
 * \param data Data struct to store info.
 */
EXTERNC ERROR_CODE recal_load_recal_info(const char *path, recal_info_t *data);

#endif /* BAM_RECAL_LIBRARY_H_ */
