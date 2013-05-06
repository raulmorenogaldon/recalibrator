#ifndef BAMAUX_H
#define BAMAUX_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "bioformats/bam-sam/alignment.h"
#include "bioformats/bam-sam/bam_file.h"

#ifdef __MMX__
#include <mmintrin.h>
#endif

#ifdef __SSE__
#include <xmmintrin.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif

//Types for sse2
//typedef char v16c __attribute__ ((vector_size (128)));

#define INFO_BATCH_SIZE 1000000000

typedef struct bam_info {
	unsigned int min_qual;
	unsigned int max_qual;
	unsigned int num_cycles;
} bam_info_t;

//Struct creation
bam_info_t *bam_new_info();

//BAM analysis
void bam_get_info(bam_file_t *bam, bam_info_t *info);
void bam_get_info_from_alignment(bam1_t* alig, bam_info_t* info);

//Struct destroy
void bam_destroy_info(bam_info_t *info);

//BAM file operations
bam_header_t* create_empty_bam_header(int num_chroms);

//Math operations
double Qsolexa(double p);
double Psolexa(double Q);
double Qsanger(double p);
double Psanger(double Q);

void compare_bams_qual(const char* bamPath0, const char* bamPath1, int cycles);

//Optimization routines
void printf_proc_features();

void print_binary(unsigned int num);

#endif
