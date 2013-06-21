#ifndef RECAL_STRUCTS_H_
#define RECAL_STRUCTS_H_

#include "config.h"

/**
 * Dinucleotide enumeration.
 * Represents all possible dinucleotide cases.
 */
enum DINUC
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
 *
 * DATA MANAGEMENT
 *
 */

/**
 * Allocate new recalibration data.
 */
recal_info_t *recal_new_info(int cycles);

/**
 * Free recalibration data.
 */
void recal_destroy_info(recal_info_t *data);

/**
 * Add recalibration data from one base.
 */
void recal_add_base(recal_info_t *data, int qual, int cycle, int dinuc, int miss);

/**
 * Add recalibration data from vector of bases
 */
void recal_add_base_v(recal_info_t *data, char *seq, char *quals, int init_cycle, int end_cycle, char *dinuc, char *misses);

/**
 * Compute deltas from bases and misses.
 */
void recal_calc_deltas(recal_info_t* data);

/**
 *
 * ENUMERATION FUNCTIONS
 *
 */

/**
 * Return dinucleotide enumeration from two bases.
 */
enum DINUC recal_get_dinuc(char A, char B);

/**
 *
 * FILE OPERATIONS
 *
 */

/**
 * Print to file data from recalibration.
 */
void recal_fprint_info(recal_info_t *data, const char *path);

/**
 * Save to file recalibration data.
 */
void recal_save_recal_info(recal_info_t *data, const char *path);

/**
 * Load from file recalibration data.
 */
void recal_load_recal_info(const char *path, recal_info_t *data);

#endif /* RECAL_STRUCTS_H_ */
