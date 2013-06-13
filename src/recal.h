#ifndef RECAL_H
#define RECAL_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "genome.h"
#include "bioformats/bam-sam/alignment.h"
#include "bioformats/bam-sam/bam_file.h"
#include "bamaux.h"
#include "timestats.h"
#include "bampool.h"

#define MIN_QUALITY 0
#define MAX_QUALITY 50

#define NUM_DINUC	17

#define SMOOTH_CONSTANT_MISS 1 //5
#define SMOOTH_CONSTANT_BASES 2 //16

#define NOT_COUNT_NUCLEOTIDE_N

#define MAX_BATCH_SIZE 100000000
//#define MAX_BATCH_SIZE 10000000

//#define USE_BATCH_POOL

//QUALITY MEASURE METHOD
//#define P_SOLEXA	//Comment to use sanger method instead solexa

//TIME MEASURES
#define D_TIME_DEBUG	//Time measurements
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
	
	
//DEBUG
#ifdef DEBUG
	//#define D_SAMPLES_OUTPUT	//Print samples to file
	
	//#define D_MAX_READS 400000	//Max reads to process (comment to ommit)
	
	//#define D_MAX_READS_W 100000	//Max reads to recal (comment to ommit)
	
	//#define D_RECAL_INTER_RESULTS
	#define D_RECAL_INTER_PROB 3300
	
	#ifdef D_SAMPLES_OUTPUT
		#define D_PROB_SAMPLE 100000   // 1/X in output samples
		char *out_path;	//Path for sample alignments file
	#endif
#endif

typedef struct recal_info {
	unsigned int min_qual;
	unsigned int num_quals;
	unsigned int num_cycles;
	unsigned int num_dinuc;
	
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

//Dinuc list
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
	d_X = 16
} dinuc;

//Unmapped alignments counter
int unmapped;

//Threads
#ifdef PTHREADS_ENABLE
typedef struct th_params {
	bam_batch_t* batch;
	genome_t* ref;
	recal_info_t* output_data;
} th_params_t;
#endif


//Allocate functions
recal_info_t *recal_new_info(int cycles);

void initialize_vector(unsigned int *vector, int size, int value);
unsigned int *new_vector(int size, int value);
double *new_vector_d(int size, double value);

//Manage functions
void recal_add_base(recal_info_t *data, int qual, int cycle, int dinuc, int miss);
void recal_add_base_v(recal_info_t *data, char *seq, char *quals, int init_cycle, int end_cycle, char *dinuc, char *misses);

//Preprocess fasta reference
void recal_preprocess_fasta(char* input_fasta, char* output_file);

//Preprocess bam data functions
void recal_get_data_from_file(char *bam_path, char *ref_name, char *ref_path, recal_info_t *out_info);
void recal_get_data_from_bam(bam_file_t *bam, genome_t* ref, recal_info_t* output_data);
void recal_get_data_from_bam_batch(bam_batch_t* batch, genome_t* ref, recal_info_t* output_data);
void recal_get_data_from_bam_alignment(bam1_t* alig, genome_t* ref, recal_info_t* output_data);
void recal_calc_deltas(recal_info_t* data);

//Recalibrate bam functions
void recal_recalibrate_bam_file(char *orig_bam_path, recal_info_t *bam_info, char *recal_bam_path);
void recal_recalibrate_bam(bam_file_t *orig_bam_f, recal_info_t *bam_info, bam_file_t *recal_bam_f);
void recal_recalibrate_batch(bam_batch_t* batch, recal_info_t *bam_info, bam_file_t *recal_bam_f);
void recal_recalibrate_alignment(bam1_t* alig, recal_info_t *bam_info, bam_file_t *recal_bam_f);

//Get dinuc
enum DINUC recal_get_dinuc(char A, char B);

//Q and P
double Qvalue(double P);
double Pvalue(double Q);

//Free functions
void recal_destroy_info(recal_info_t *data);

//File print info
void recal_fprint_info(recal_info_t *data, const char *path);

//File recal_info save
void recal_save_recal_info(recal_info_t *data, const char *path);
void recal_load_recal_info(const char *path, recal_info_t *data);

#endif
