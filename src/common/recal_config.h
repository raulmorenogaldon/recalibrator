/*
 * recal_config.h
 *
 *  Created on: Jun 25, 2013
 *      Author: rmoreno
 */

#ifndef RECAL_CONFIG_H_
#define RECAL_CONFIG_H_

/**
 * Recalibration data struct defines
 */
#define MIN_QUALITY 0
#define MAX_QUALITY 50
#define NUM_DINUC	17

//Smoothing constants
#define SMOOTH_CONSTANT_MISS 1 //5
#define SMOOTH_CONSTANT_BASES 1 //16

/**
 * Recalibration parameters
 */
#define MIN_QUALITY_TO_STAT 6
#define NOT_COUNT_NUCLEOTIDE_N
#define STAT_FIRST_DINUC
#define NOT_MAPPING_QUAL_ZERO //Reads with map score of 0 wont be stated

/**
 * BAM management
 */
#define MAX_BATCH_SIZE 1000000000

/**
 * Time measures
 */
#define D_TIME_DEBUG
#ifdef D_TIME_DEBUG
	#define D_SLOT_GET_DATA_BAM 0
	#define D_SLOT_GET_DATA_BATCH 1
	#define D_SLOT_GET_DATA_ALIG 2
	#define D_SLOT_READ_BATCH 3
	#define D_SLOT_CALC_DELTAS 4
	#define D_SLOT_RECALIBRATE 5
	#define D_SLOT_MEMCOPY_BATCH 6
	#define D_SLOT_RECAL_ALIG 7
#endif

#endif /* RECAL_CONFIG_H_ */
