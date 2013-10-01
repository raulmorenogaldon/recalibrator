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
#define MAX_QUALITY 94
#define NUM_DINUC	17

//Smoothing constants
//#define SMOOTH_CONSTANT_MISS 1 //5
//#define SMOOTH_CONSTANT_BASES 1 //16
#define SMOOTH_CONSTANT 1

/**
 * Recalibration parameters
 */
#define MIN_QUALITY_TO_STAT 6
#define NOT_COUNT_NUCLEOTIDE_N
#define STAT_FIRST_DINUC
#define NOT_MAPPING_QUAL_ZERO //Reads with map score of 0 wont be stated
#define NOT_PRIMARY_ALIGNMENT
//#define CHECK_DUPLICATES
//#define SPLIT_BATCHS_BY_CHROM
#define USE_SSE

//#define D_MAX_READS_W 1000

/**
 * BAM management
 */
#define MAX_BATCH_SIZE 1000000000
//#define MAX_BATCH_SIZE 1000000

/**
 * Probability math limits
 */
#define P_SOLEXA_MAX 62
#define P_SOLEXA_MIN -5
#define P_SANGER_MAX 93
#define P_SANGER_MIN 0

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
	#define D_SLOT_WRITE_BATCH 8
#endif

#endif /* RECAL_CONFIG_H_ */
