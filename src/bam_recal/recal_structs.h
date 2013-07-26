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
	qual_t min_qual;	//Set minor stored quality
	uint32_t num_quals;	//Set range of stored qualities
	uint32_t num_cycles;	//Set maximum number of cycles stored
	uint32_t num_dinuc;	//Set number of maximum dinucleotides

	error_t total_miss;				//Total misses
	uint32_t total_bases;				//Total bases
	delta_t total_delta;			//Global delta

    error_t* qual_miss;				//Misses per quality
    uint32_t* qual_bases;				//Bases per quality
    delta_t* qual_delta;			//Delta per quality

    error_t* qual_cycle_miss;		//Misses per quality-cycle pair
    uint32_t* qual_cycle_bases;		//Bases per quality-cycle pair
    delta_t* qual_cycle_delta;		//Deltas per quality-cycle pair

    error_t* qual_dinuc_miss;		//Misses per quality-dinuc pair
    uint32_t* qual_dinuc_bases;		//Bases per quality-dinuc pair
    delta_t* qual_dinuc_delta;		//Deltas per quality-dinuc pair
};

/**
 * PRIVATE FUNCTIONS
 */

#endif /* RECAL_STRUCTS_H_ */
