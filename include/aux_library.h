/*
 * aux_library.h
 *
 *  Created on: Jun 28, 2013
 *      Author: rmoreno
 */

#ifndef AUX_LIBRARY_H_
#define AUX_LIBRARY_H_

//#include "bam_recal_library.h"
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

EXTERNC double gaussian_function(double value, double a, double b, double c, double d);

EXTERNC double log10_gamma(uint32_t n);

/***************************
 * QUALITY OPERATIONS
 **************************/

EXTERNC ERROR_CODE recal_get_estimated_Q(uint32_t *v_bases, uint32_t count, uint8_t start_quality, double *estimated_Q);

EXTERNC ERROR_CODE recal_get_empirical_Q(double miss, uint32_t bases, double initial_quality, double *emp_qual);

EXTERNC ERROR_CODE log10_Qemp_Reported(double Qemp, double Qreported, double *log);

EXTERNC ERROR_CODE log10_Qemp_likelihood(double Qemp, uint32_t obs, uint32_t err, double *log);

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

EXTERNC char * new_sequence_from_bam(bam1_t *bam1);

EXTERNC char * new_quality_from_bam(bam1_t *bam1, int base_quality);

EXTERNC ERROR_CODE decompose_cigar(char *cigar, uint8_t cigar_length, char *n_elem, char *type, uint8_t *types_length, uint8_t max_types_length);

/***************************
 * VECTOR OPERATIONS
 **************************/

/**
 * Initializes vector with initial values.
 * \param vector Vector to initialize
 * \param size Size of vector
 * \param value Value to initialize vector
 */
EXTERNC ERROR_CODE initialize_vector(uint32_t *vector, const size_t size, const uint32_t value);

/**
 * Return vector of integers with initial values.
 * \param size Size of the new vector.
 * \param value Initial value in all elements of the vector
 * \param out_vector Pointer to pointer which stores the vector
 */

/**
 * Return vector of double with initial values.
 * \param size Size of the new vector.
 * \param value Initial value in all elements of the vector
 * \param out_vector Pointer to pointer which stores the vector
 */
EXTERNC ERROR_CODE new_vector_uint32(const size_t size, const uint32_t value, uint32_t **out_vector);

EXTERNC ERROR_CODE new_vector_double(const size_t size, const double value, double **out_vector);

EXTERNC ERROR_CODE max_value(double *vector, size_t size, double *max);

EXTERNC ERROR_CODE max_index(double *vector, size_t size, uint16_t *max_i);

/***************************
 * MISCELANEA OPERATIONS
 **************************/

EXTERNC void printf_proc_features();

EXTERNC void print_binary(unsigned int num);

#endif /* AUX_LIBRARY_H_ */
