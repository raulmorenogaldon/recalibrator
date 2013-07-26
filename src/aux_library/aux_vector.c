#include "aux_vector.h"

/**
 * Initializes vector with initial values.
 */
ERROR_CODE
initialize_vector(uint32_t *vector, const size_t size, const uint32_t value)
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
new_vector_bases(const size_t size, const base_t value, base_t **out_vector)
{
	int i;
	base_t *vector;

	if(size == 0)
	{
		*out_vector = NULL;
		return INVALID_INPUT_SIZE_0;
	}

	vector = (base_t *)malloc(size * sizeof(base_t));

	for(i = 0; i < size; i++)
		vector[i] = value;

	printf("Created new vector with %d positions, total size %lu bytes\n", size, size * sizeof(base_t));

	*out_vector = vector;

	return NO_ERROR;
}

/**
 * Return vector of double with initial values.
 */
ERROR_CODE
new_vector_miss(const size_t size, const error_t value, error_t **out_vector)
{
	int i;
	error_t *vector;

	if(size == 0)
	{
		*out_vector = NULL;
		return INVALID_INPUT_SIZE_0;
	}

	vector = (error_t *)malloc(size * sizeof(error_t));

	for(i = 0; i < size; i++)
		vector[i] = value;

	printf("Created new vector with %d positions, total size %lu bytes\n", size, size * sizeof(error_t));

	*out_vector = vector;

	return NO_ERROR;
}

ERROR_CODE
new_vector_delta(const size_t size, const delta_t value, delta_t **out_vector)
{
	int i;
	delta_t *vector;

	if(size == 0)
	{
		*out_vector = NULL;
		return INVALID_INPUT_SIZE_0;
	}

	vector = (delta_t *)malloc(size * sizeof(delta_t));

	for(i = 0; i < size; i++)
		vector[i] = value;

	printf("Created new vector with %d positions, total size %lu bytes\n", size, size * sizeof(delta_t));

	*out_vector = vector;

	return NO_ERROR;
}

ERROR_CODE
max_value(double *vector, size_t size, double *max)
{
	double maximum = -DBL_MAX;
	double aux;
	int i;

	for(i = size-1; i >= 0; i--)
	{
		aux = vector[i];
		if(aux > maximum)
			maximum = aux;
	}

	*max = maximum;

	return NO_ERROR;
}

ERROR_CODE
max_index(double *vector, size_t size, uint16_t *max_i)
{
	double max = -DBL_MAX;
	int index = -1;
	double aux;
	int i;

	for(i = size-1; i >= 0; i--)
	{
		aux = vector[i];
		if(aux > max)
		{
			max = aux;
			index = i;
		}
	}

	*max_i = index;

	return NO_ERROR;
}
