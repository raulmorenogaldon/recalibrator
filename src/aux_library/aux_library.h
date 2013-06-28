/*
 * aux_library.h
 *
 *  Created on: Jun 28, 2013
 *      Author: rmoreno
 */

#ifndef AUX_LIBRARY_H_
#define AUX_LIBRARY_H_

#include "bam.h"


/***************************
 * MATH OPERATIONS
 **************************/

/**
 * Returns quality from probability.
 */
extern inline double Qvalue(double P);

/**
 * Return probability from quality.
 */
extern inline double Pvalue(double Q);

/**
 * Return Solexa quality from probability.
 */
extern inline double Qsolexa(double p);

/**
 * Return probability from Solexa quality.
 */
extern inline double Psolexa(double Q);

/**
 * Return Sanger quality from probability.
 */
extern inline double Qsanger(double p);

/**
 * Return probability from Solexa quality.
 */
extern inline double Psanger(double Q);


/***************************
 * BAM OPERATIONS
 **************************/

extern bam_header_t* create_empty_bam_header(int num_chroms);

extern void compare_bams_qual(const char* bamPath0, const char* bamPath1, int cycles);


/***************************
 * VECTOR OPERATIONS
 **************************/

/**
 * Initializes vector with initial values.
 */
extern void initialize_vector(unsigned int *vector, int size, int value);

/**
 * Return vector of integers with initial values.
 */
extern unsigned int *new_vector(int size, int value);

/**
 * Return vector of double with initial values.
 */
extern double *new_vector_d(int size, double value);


/***************************
 * MISCELANEA OPERATIONS
 **************************/

extern void printf_proc_features();

extern void print_binary(unsigned int num);

#endif /* AUX_LIBRARY_H_ */
