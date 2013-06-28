/*
 * bam_recal_library.h
 *
 *  Created on: Jun 28, 2013
 *      Author: rmoreno
 */

#ifndef BAM_RECAL_LIBRARY_H_
#define BAM_RECAL_LIBRARY_H_

#include "bam.h"
#include "bam_file.h"
#include "alignment.h"
#include "genome.h"

/**
 * Recalibration data storage
 * This struct hold all data necessary to perform recalibrations
 */
typedef struct recal_info {
	unsigned int min_qual;	//Set minor stored quality
	unsigned int num_quals;	//Set range of stored qualities
	unsigned int num_cycles;	//Set maximum number of cycles stored
	unsigned int num_dinuc;	//Set number of maximum dinucleotides

    long unsigned int total_miss;	//Total misses
    long unsigned int total_bases;	//Total bases
    double total_delta;	//Global delta
    unsigned int* qual_miss;	//Misses per quality
    unsigned int* qual_bases;	//Bases per quality
    double* qual_delta;	//Delta per quality
    unsigned int* qual_cycle_miss;		//Misses per quality-cycle pair
    unsigned int* qual_cycle_bases;		//Bases per quality-cycle pair
    double* qual_cycle_delta;		//Deltas per quality-cycle pair
    unsigned int* qual_dinuc_miss;		//Misses per quality-dinuc pair
    unsigned int* qual_dinuc_bases;		//Bases per quality-dinuc pair
    double* qual_dinuc_delta;		//Deltas per quality-dinuc pair
} recal_info_t;

/**
 * Dinucleotide enumeration.
 * Represents all possible dinucleotide cases.
 */
extern enum DINUC
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
};


/***********************************************
 * DATA MANAGEMENT
 **********************************************/

/**
 * Allocate new recalibration data.
 */
extern recal_info_t *recal_new_info(int cycles);

/**
 * Free recalibration data.
 */
extern void recal_destroy_info(recal_info_t *data);

/**
 * Add recalibration data from one base.
 */
extern void recal_add_base(recal_info_t *data, int qual, int cycle, int dinuc, int miss);

/**
 * Add recalibration data from vector of bases
 */
extern void recal_add_base_v(recal_info_t *data, char *seq, char *quals, int init_cycle, int end_cycle, char *dinuc, char *misses);

/**
 * Compute deltas from bases and misses.
 */
extern void recal_calc_deltas(recal_info_t* data);


/***********************************************
 * DINUC ENUM OPERATIONS
 **********************************************/

/**
 * Return dinucleotide enumeration from two bases.
 */
extern enum DINUC recal_get_dinuc(char A, char B);


/***********************************************
 * BAM RECALIBRATION PHASE 1 - DATA COLLECT
 **********************************************/

/**
 * Get recalibration data from BAM path.
 */
extern void recal_get_data_from_file(char *bam_path, char *ref_name, char *ref_path, recal_info_t *out_info);

/**
 * Get recalibration data from BAM file.
 */
extern void recal_get_data_from_bam(bam_file_t *bam, genome_t* ref, recal_info_t* output_data);

/**
 * Get recalibration data from BAM batch of alignments.
 */
extern void recal_get_data_from_bam_batch(bam_batch_t* batch, genome_t* ref, recal_info_t* output_data);

/**
 * Get recalibration data from alignment.
 */
extern void recal_get_data_from_bam_alignment(bam1_t* alig, genome_t* ref, recal_info_t* output_data);


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


/***********************************************
 * FILE OPERATIONS
 **********************************************/

/**
 * Print to file data from recalibration.
 */
extern void recal_fprint_info(recal_info_t *data, const char *path);

/**
 * Save to file recalibration data.
 */
extern void recal_save_recal_info(recal_info_t *data, const char *path);

/**
 * Load from file recalibration data.
 */
extern void recal_load_recal_info(const char *path, recal_info_t *data);

#endif /* BAM_RECAL_LIBRARY_H_ */
