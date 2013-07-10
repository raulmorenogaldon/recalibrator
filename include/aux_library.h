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
EXTERNC double Qvalue(double P);

/**
 * Return probability from quality.
 */
EXTERNC double Pvalue(double Q);

/**
 * Return Solexa quality from probability.
 */
EXTERNC double Qsolexa(double p);

/**
 * Return probability from Solexa quality.
 */
EXTERNC double Psolexa(double Q);

/**
 * Return Sanger quality from probability.
 */
EXTERNC double Qsanger(double p);

/**
 * Return probability from Solexa quality.
 */
EXTERNC double Psanger(double Q);


/***************************
 * BAM OPERATIONS
 **************************/

EXTERNC ERROR_CODE create_empty_bam_header(const unsigned int num_chroms, bam_header_t *header);

EXTERNC ERROR_CODE compare_bams_qual(const char* bamPath0, const char* bamPath1, const int cycles);


/***************************
 * VECTOR OPERATIONS
 **************************/

/**
 * Initializes vector with initial values.
 */
EXTERNC ERROR_CODE initialize_vector(unsigned int *vector, const size_t size, const int value);

/**
 * Return vector of integers with initial values.
 */
EXTERNC ERROR_CODE new_vector(const size_t size, const int value, unsigned int *out_vector);

/**
 * Return vector of double with initial values.
 */
EXTERNC ERROR_CODE new_vector_d(const size_t size, const double value, double *out_vector);


/***************************
 * MISCELANEA OPERATIONS
 **************************/

EXTERNC void printf_proc_features();

EXTERNC void print_binary(unsigned int num);

#endif /* AUX_LIBRARY_H_ */
