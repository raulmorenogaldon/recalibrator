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
 * \param p Error probability [0,1]
 */
EXTERNC double Qvalue(double p);

/**
 * Return probability from quality.
 * \param q Error quality. Default -> [P_SANGER_MIN,P_SANGER_MAX]
 */
EXTERNC double Pvalue(double q);

/**
 * Return Solexa quality from probability.
 * \param p Error probability [0,1]
 */
EXTERNC double Qsolexa(double p);

/**
 * Return probability from Solexa quality.
 * \param q Error quality [P_SOLEXA_MIN,P_SOLEXA_MAX]
 */
EXTERNC double Psolexa(double q);

/**
 * Return Sanger quality from probability.
 * \param p Error probability [0,1]
 */
EXTERNC double Qsanger(double p);

/**
 * Return probability from Sanger quality.
 * \param q Error quality [P_SANGER_MIN,P_SANGER_MAX]
 */
EXTERNC double Psanger(double q);


/***************************
 * BAM OPERATIONS
 **************************/

/**
 * Init bam_header_t struct to work with recalibration.
 * \param num_chroms Number of target chromosomes
 * \param header Pointer to previously allocated empty header struct
 */
EXTERNC ERROR_CODE init_empty_bam_header(const unsigned int num_chroms, bam_header_t *header);

EXTERNC ERROR_CODE compare_bams_qual(const char* bamPath0, const char* bamPath1, const int cycles);


/***************************
 * VECTOR OPERATIONS
 **************************/

/**
 * Initializes vector with initial values.
 * \param vector Vector to initialize
 * \param size Size of vector
 * \param value Value to initialize vector
 */
EXTERNC ERROR_CODE initialize_vector(unsigned int *vector, const size_t size, const unsigned int value);

/**
 * Return vector of integers with initial values.
 * \param size Size of the new vector.
 * \param value Initial value in all elements of the vector
 * \param out_vector Pointer to pointer which stores the vector
 */
EXTERNC ERROR_CODE new_vector(const size_t size, const int value, unsigned int **out_vector);

/**
 * Return vector of double with initial values.
 * \param size Size of the new vector.
 * \param value Initial value in all elements of the vector
 * \param out_vector Pointer to pointer which stores the vector
 */
EXTERNC ERROR_CODE new_vector_d(const size_t size, const double value, double **out_vector);


/***************************
 * MISCELANEA OPERATIONS
 **************************/

EXTERNC void printf_proc_features();

EXTERNC void print_binary(unsigned int num);

#endif /* AUX_LIBRARY_H_ */
