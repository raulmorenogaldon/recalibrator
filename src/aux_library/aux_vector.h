#ifndef AUX_VECTOR_H_
#define AUX_VECTOR_H_

#include <stdio.h>
#include <stdlib.h>
#include <recal_common.h>

/**
 * Initializes vector with initial values.
 */
ERROR_CODE initialize_vector(unsigned int *vector, const size_t size, const int value);

/**
 * Return vector of integers with initial values.
 */
ERROR_CODE new_vector(const size_t size, const int value, unsigned int *out_vector);

/**
 * Return vector of double with initial values.
 */
ERROR_CODE new_vector_d(const size_t size, const double value, double *out_vector);

#endif /* AUX_VECTOR_H_ */
