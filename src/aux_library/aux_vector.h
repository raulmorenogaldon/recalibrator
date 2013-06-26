#ifndef AUX_VECTOR_H_
#define AUX_VECTOR_H_

#include <stdio.h>
#include <stdlib.h>

#include "common.h"

/**
 * Initializes vector with initial values.
 */
void initialize_vector(unsigned int *vector, int size, int value);

/**
 * Return vector of integers with initial values.
 */
unsigned int *new_vector(int size, int value);

/**
 * Return vector of double with initial values.
 */
double *new_vector_d(int size, double value);

#endif /* AUX_VECTOR_H_ */
