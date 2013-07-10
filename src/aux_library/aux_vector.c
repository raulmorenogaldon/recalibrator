#include "aux_vector.h"

/**
 * Initializes vector with initial values.
 */
ERROR_CODE
initialize_vector(unsigned int *vector, const size_t size, const int value)
{
	int i;

	if(size == 0)
		return INVALID_INPUT_SIZE_0;

	for(i = 0; i < size; i++)
	{
		vector[i] = value;
	}

	return NO_ERROR;
}

/**
 * Return vector of integers with initial values.
 */
ERROR_CODE
new_vector(const size_t size, const int value, unsigned int *out_vector)
{
	int i;

	if(size == 0)
		return INVALID_INPUT_SIZE_0;

	out_vector = (unsigned int *)malloc(size * sizeof(unsigned int));

	for(i = 0; i < size; i++)
		out_vector[i] = value;

	printf("Created new vector with %d positions, total size %lu bytes\n", size, size * sizeof(unsigned int));

	return NO_ERROR;
}

/**
 * Return vector of double with initial values.
 */
ERROR_CODE
new_vector_d(const size_t size, const double value, double *out_vector)
{
	int i;

	if(size == 0)
		return INVALID_INPUT_SIZE_0;

	out_vector = (double *)malloc(size * sizeof(double));

	for(i = 0; i < size; i++)
		out_vector[i] = value;

	printf("Created new vector with %d positions, total size %lu bytes\n", size, size * sizeof(double));

	return NO_ERROR;
}
