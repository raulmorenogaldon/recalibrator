#include "aux_vector.h"

/**
 * Initializes vector with initial values.
 */
void initialize_vector(unsigned int *vector, int size, int value)
{
	int i;

	for(i = 0; i < size; i++)
	{
		vector[i] = value;
	}
}

/**
 * Return vector of integers with initial values.
 */
unsigned int *new_vector(unsigned int size, int value)
{
	unsigned int *vector;
	int i;

	vector = (unsigned int *)malloc(size * sizeof(unsigned int));

	for(i = 0; i < size; i++)
		vector[i] = value;

	printf("Created new vector with %d positions, total size %lu bytes\n", size, size * sizeof(unsigned int));

	return vector;
}

/**
 * Return vector of double with initial values.
 */
double *new_vector_d(unsigned int size, double value)
{
	double *vector;
	int i;

	vector = (double *)malloc(size * sizeof(double));

	for(i = 0; i < size; i++)
		vector[i] = value;

	printf("Created new vector with %d positions, total size %lu bytes\n", size, size * sizeof(double));

	return vector;
}
