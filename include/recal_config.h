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
//#define MAX_BATCH_SIZE 	1000000000
#define MAX_BATCH_SIZE 	100000000
//#define MAX_BATCH_SIZE 	10000000
//#define MAX_BATCH_SIZE 	1000000

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
	//DELTAS
	#define D_SLOT_PROCCESS_DELTAS 0

	//PHASE 1
	#define D_SLOT_PH1_COLLECT_BAM 1
	#define D_SLOT_PH1_READ_BATCH 2
	#define D_SLOT_PH1_COLLECT_BATCH 3
	#define D_SLOT_PH1_COLLECT_ALIG 4
	#define D_SLOT_PH1_COLLECT_REDUCE_DATA 10
	#define D_SLOT_PH1_ITERATION 11

	//PHASE 2
	#define D_SLOT_PH2_READ_BATCH 5
	#define D_SLOT_PH2_PROCCESS_BATCH 6
	#define D_SLOT_PH2_WRITE_BATCH 7
	#define D_SLOT_PH2_RECAL_ALIG 8
	#define D_SLOT_PH2_RECALIBRATE 9
	#define D_SLOT_PH2_ITERATION 12
#endif

/**
 * OpenMP
 */
//#define D_TIME_OPENMP_VERBOSE
//#define OMP_SCHEDULE "dynamic"

#endif /* RECAL_CONFIG_H_ */
