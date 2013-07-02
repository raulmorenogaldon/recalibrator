#ifndef AUX_VECTOR_H_
#define AUX_VECTOR_H_

#include <stdio.h>
#include <stdlib.h>
#include <recal_common.h>

/**
 * Initializes vector with initial values.
 */
extern void initialize_vector(unsigned int *vector, unsigned int size, int value);

/**
 * Return vector of integers with initial values.
 */
extern unsigned int *new_vector(unsigned int size, int value);

/**
 * Return vector of double with initial values.
 */
extern double *new_vector_d(unsigned int size, double value);

#endif /* AUX_VECTOR_H_ */
