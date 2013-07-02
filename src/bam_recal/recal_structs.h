#ifndef RECAL_STRUCTS_H_
#define RECAL_STRUCTS_H_

#include <recal_common.h>
#include <recal_config.h>
#include <aux_library.h>
#include <timestats.h>

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

/**
 * Return dinucleotide enumeration from two bases.
 */
extern enum DINUC recal_get_dinuc(char A, char B);


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

#endif /* RECAL_STRUCTS_H_ */
