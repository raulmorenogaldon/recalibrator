#ifndef RECAL_STRUCTS_H_
#define RECAL_STRUCTS_H_

#include <recal_common.h>
#include <recal_config.h>
#include <aux_library.h>
#include <bam_recal_library.h>
#include <timestats.h>

/**
 * Recalibration data storage
 * This struct hold all data necessary to perform recalibrations
 */
struct recal_info {
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
};

/**
 * PRIVATE FUNCTIONS
 */

#endif /* RECAL_STRUCTS_H_ */
