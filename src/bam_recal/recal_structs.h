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
	uint8_t min_qual;	//Set minor stored quality
	uint32_t num_quals;	//Set range of stored qualities
	uint32_t num_cycles;	//Set maximum number of cycles stored
	uint32_t num_dinuc;	//Set number of maximum dinucleotides

	double total_miss;				//Total misses
	uint32_t total_bases;				//Total bases
	double total_delta;			//Global delta

	double* qual_miss;				//Misses per quality
    uint32_t* qual_bases;				//Bases per quality
    double* qual_delta;			//Delta per quality

    double* qual_cycle_miss;		//Misses per quality-cycle pair
    uint32_t* qual_cycle_bases;		//Bases per quality-cycle pair
    double* qual_cycle_delta;		//Deltas per quality-cycle pair

    double* qual_dinuc_miss;		//Misses per quality-dinuc pair
    uint32_t* qual_dinuc_bases;		//Bases per quality-dinuc pair
    double* qual_dinuc_delta;		//Deltas per quality-dinuc pair
};


typedef struct data_collect_env {

	//Sequence storage
	char *bam_seq;
	char *bam_quals;

	//Auxiliar
	char *aux_res_seq;
	char *aux_res_qual;

	//Maximum length
	uint32_t bam_seq_max_l;

} recal_data_collect_env_t;

/**
 * PRIVATE FUNCTIONS
 */

#endif /* RECAL_STRUCTS_H_ */
