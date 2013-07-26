#include "recal_structs.h"

#include <math.h>
#include <float.h>
#include <fenv.h>

/**
 *
 * DATA MANAGEMENT
 *
 */

/**
 * Allocate new recalibration data.
 */
ERROR_CODE
recal_init_info(const uint32_t cycles, recal_info_t **out_data)
{
	recal_info_t *data;
	int vector_size = MAX_QUALITY - MIN_QUALITY;

	data = (recal_info_t *)malloc(sizeof(recal_info_t));

	//Struct initialization
	data->min_qual = MIN_QUALITY;
	data->num_quals = vector_size;
	data->num_cycles = cycles;
	data->num_dinuc = NUM_DINUC;

	//Total counters
	data->total_miss = 0.0;
	data->total_bases = 0;
	data->total_delta = 0.0;

	//Quality vectors
	new_vector_miss(vector_size, 0.0, &(data->qual_miss));
	new_vector_bases(vector_size, 0.0, &(data->qual_bases));
	new_vector_delta(vector_size, 0.0, &(data->qual_delta));

	//Qual-Cycle matrix
	new_vector_miss(vector_size * cycles, 0.0, &(data->qual_cycle_miss));
	new_vector_bases(vector_size * cycles, 0.0, &(data->qual_cycle_bases));
	new_vector_delta(vector_size * cycles, 0.0, &(data->qual_cycle_delta));

	//Qual-Dinuc matrix
	new_vector_miss(vector_size * cycles, 0.0, &(data->qual_dinuc_miss));
	new_vector_bases(vector_size * cycles, 0.0, &(data->qual_dinuc_bases));
	new_vector_delta(vector_size * cycles, 0.0, &(data->qual_dinuc_delta));

	*out_data = data;

	return NO_ERROR;
}

/**
 * Free recalibration data.
 */
ERROR_CODE
recal_destroy_info(recal_info_t **data)
{
	recal_info_t *d = *data;

	//Free quality vector
	free(d->qual_miss);
	free(d->qual_bases);
	free(d->qual_delta);

	//Free quality-cycle matrix
	free(d->qual_cycle_miss);
	free(d->qual_cycle_bases);
	free(d->qual_cycle_delta);

	//Free quality-dinuc matrix
	free(d->qual_dinuc_miss);
	free(d->qual_dinuc_bases);
	free(d->qual_dinuc_delta);

	//Free struct
	free(d);
	*data = NULL;

	return NO_ERROR;
}

/**
 * Add recalibration data from one base.
 */
ERROR_CODE
recal_add_base(recal_info_t *data, const qual_t qual, const cycle_t cycle, const dinuc_t dinuc, const error_t match)
{
	int qual_index = qual - data->min_qual;
	int qual_cycle_index = qual_index * data->num_cycles + cycle;
	int qual_dinuc_index = qual_index * data->num_dinuc + dinuc;

	if(qual < MIN_QUALITY_TO_STAT)
		return INVALID_INPUT_QUAL;

	//Error check
	if(qual_index >= MAX_QUALITY - MIN_QUALITY || qual_index < 0)
	{
		printf("add_base: ERROR, qual must be positive and minor than MAX_QUALITY - MIN_QUALITY ==> Qual = %d\n", qual);
		return INVALID_INPUT_QUAL;
	}
	if(cycle < 0 || cycle > data->num_cycles - 1)
	{
		printf("add_base: ERROR, cycle must be positive and minor than NUM_CYCLES\n");
		return INVALID_INPUT_QUAL;
	}
	if(dinuc > NUM_DINUC - 1)
	{
		printf("add_base: ERROR, dinuc must be minor than NUM_DINUC %d \n", dinuc);
		return INVALID_INPUT_QUAL;
	}

	//Increase total counters
	if(!match)
	data->total_miss++;
	data->total_bases++;

	//Increase quality vector
	if(!match)
	data->qual_miss[qual_index]++;
	data->qual_bases[qual_index]++;

	//Increase quality-cycle matrix
	if(!match)
	data->qual_cycle_miss[qual_cycle_index]++;
	data->qual_cycle_bases[qual_cycle_index]++;

	//Increase quality-dinuc matrix (if not dinuc != -1)
	if(dinuc >= 0)
	{
		if(!match)
		data->qual_dinuc_miss[qual_dinuc_index]++;
		data->qual_dinuc_bases[qual_dinuc_index]++;
	}
	else
	{
		printf("ERROR: unrecognized dinuc Q: %d, C: %d, D: %d, Match: %d\n", qual, cycle, dinuc, match);
		return INVALID_DINUCLEOTIDE;
	}

	return NO_ERROR;
}

/**
 * Add recalibration data from vector of bases
 */
ERROR_CODE
recal_add_base_v(recal_info_t *data, const base_t *seq, const qual_t *quals, const cycle_t init_cycle, const uint32_t num_cycles, const dinuc_t *dinuc, const error_t *matches)
{
	int i;
	uint8_t qual_index;
	uint16_t qual_cycle_index;
	uint16_t qual_dinuc_index;
	uint32_t cycles;

	cycles = num_cycles;

	//Iterates cycles
	//for(i = init_cycle; i <= end_cycle; i++)
	for(i = 0; i < cycles; i++)
	{
		switch(seq[i])
		{
		case 'A':
		case 'C':
		case 'G':
		case 'T':
		#ifndef NOT_COUNT_NUCLEOTIDE_N
		case 'N':
		#endif
			recal_add_base(data, quals[i], i + init_cycle, dinuc[i], matches[i]);
			break;

		default:
			//printf("ERROR: Corrupted nucleotide read = %c\n", seq[i]);
			break;
		}

		/*#ifdef NOT_COUNT_NUCLEOTIDE_N
		if(seq[i] != 'N' && quals[i - init_cycle] - 1 >= MIN_QUALITY_TO_STAT)
		#endif
		{
			//Indices
			qual_index = quals[i] - data->min_qual - 1;
			qual_cycle_index = qual_index * data->num_cycles + i;
			qual_dinuc_index = qual_index * data->num_dinuc + dinuc[i];

			//Increase bases
			data->total_bases++;
			data->qual_bases[qual_index]++;
			data->qual_cycle_bases[qual_cycle_index]++;
			if(dinuc >= 0)
			{
				data->qual_dinuc_bases[qual_dinuc_index]++;
			}

			//Increase misses
			if(!matches[i])
			{
				data->total_miss++;
				data->qual_miss[qual_index]++;
				data->qual_cycle_miss[qual_cycle_index]++;
				if(dinuc >= 0)
				{
					data->qual_dinuc_miss[qual_dinuc_index]++;
				}
			}
		}*/
	}

	return NO_ERROR;
}

/**
 * Compute deltas from bases and misses.
 */
ERROR_CODE
recal_calc_deltas(recal_info_t* data)
{
	double global_empirical, phred, err0;
	//double r_empirical;
	int matrix_index;
	int i, j;

	double estimated_Q;
	double emp_Q;

	//Time measures
	#ifdef D_TIME_DEBUG
		time_init_slot(D_SLOT_CALC_DELTAS, clock(), TIME_GLOBAL_STATS);
	#endif

	printf("Processing deltas...\n");

	//Estimated Q
	recal_get_estimated_Q(data->qual_bases, data->num_quals, 0, &estimated_Q);

	//Global delta
	//data->total_miss = 11;
	//global_empirical = (double)(data->total_miss) / (double)(data->total_bases);
	//phred = 0.0;
	//sum_errors = 0.0;
	//for(i = 0; i < data->num_quals; i++)
	{
		//if(data->min_qual + i != 0)
		{
			//err0 = 1.0 / (10.0 * log10(data->min_qual + i + 1));
			//err0 = 1.0 / ((double)(data->min_qual + i));
			//err0 = Pvalue((double)(/*data->min_qual +*/ i));
			//phred += err0 * ((double)data->qual_bases[i] - 1.0);
		}
	}
	//phred = phred / (double)data->total_bases;
	//int calidad_phred =  Qvalue(phred);
	//int calidad_global = Qvalue(global_empirical);
	//estimated_Q = Qvalue(sum_errors / (double)data->total_bases);

	//data->total_delta = Qvalue(global_empirical) - Qvalue(phred);

	//Get empirical global qual
	recal_get_empirical_Q(data->total_miss, data->total_bases, estimated_Q, &emp_Q);

	//Calc global delta
	data->total_delta = emp_Q - estimated_Q;


	//Delta R
	for(i = 0; i < data->num_quals; i++)
	{
		if(data->qual_bases[i] != 0)
		{
			//int calidad = Qvalue( (double)(data->qual_miss[i]) / (double)(data->qual_bases[i]) );
			//data->qual_delta[i] = calidad
			//		- (double)(i /*+ data->min_qual*/)
			//		- data->total_delta;

			recal_get_empirical_Q(data->qual_miss[i], data->qual_bases[i], data->total_delta + estimated_Q, &emp_Q);
			data->qual_delta[i] = emp_Q - (data->total_delta + estimated_Q);
		}
	}


	//Delta R,C
	for(i = 0; i < data->num_quals; i++)
	{
		for(j = 0; j < data->num_cycles; j++)
		{
			matrix_index = i * data->num_cycles + j;
			if(data->qual_cycle_bases[matrix_index] != 0)
			{
				//data->qual_cycle_delta[matrix_index] = Qvalue((double)(data->qual_cycle_miss[matrix_index]) / (double)(data->qual_cycle_bases[matrix_index]))
				//	- (double)(i /*+ data->min_qual*/)
				//	- (data->total_delta + data->qual_delta[i]);
				recal_get_empirical_Q(data->qual_cycle_miss[matrix_index], data->qual_cycle_bases[matrix_index],
						data->qual_delta[i] + data->total_delta + estimated_Q, &emp_Q);
				data->qual_cycle_delta[i] = emp_Q - (data->qual_delta[i] + data->total_delta + estimated_Q);
			}
		}
	}

	//Delta R,D
	for(i = 0; i < data->num_quals; i++)
	{
		for(j = 0; j < data->num_dinuc; j++)
		{
			matrix_index = i * data->num_dinuc + j;
			if(data->qual_dinuc_bases[matrix_index] != 0)
			{
				//data->qual_dinuc_delta[matrix_index] = Qvalue((double)(data->qual_dinuc_miss[matrix_index]) / (double)(data->qual_dinuc_bases[matrix_index]))
				//	- (double)(i /*+ data->min_qual*/)
				//	- (data->total_delta + data->qual_delta[i]);
				recal_get_empirical_Q(data->qual_dinuc_miss[matrix_index], data->qual_dinuc_bases[matrix_index],
						data->qual_delta[i] + data->total_delta + estimated_Q, &emp_Q);
				data->qual_dinuc_delta[i] = emp_Q - (data->qual_delta[i] + data->total_delta + estimated_Q);
			}
		}
	}

	printf("Deltas processed.\n");

	#ifdef D_TIME_DEBUG
		time_set_slot(D_SLOT_CALC_DELTAS, clock(), TIME_GLOBAL_STATS);
	#endif

	return NO_ERROR;
}

/**
 *
 * ENUMERATION FUNCTIONS
 *
 */

/**
 * Return dinucleotide enumeration from two bases.
 */
ERROR_CODE
recal_get_dinuc(const char A, const char B, dinuc_t *out_dinuc)
{
	*out_dinuc = d_X;

	switch(A)
	{
			case 'A':
			switch(B)
			{
					case 'A':
						*out_dinuc = dAA;
					break;
					case 'G':
						*out_dinuc = dAG;
					break;
					case 'C':
						*out_dinuc = dAC;
					break;
					case 'T':
						*out_dinuc = dAT;
					break;
			}
			break;
			case 'G':
			switch(B)
			{
					case 'A':
						*out_dinuc = dGA;
					break;
					case 'G':
						*out_dinuc = dGG;
					break;
					case 'C':
						*out_dinuc = dGC;
					break;
					case 'T':
						*out_dinuc = dGT;
					break;
			}
			break;
			case 'C':
			switch(B)
			{
					case 'A':
						*out_dinuc = dCA;
					break;
					case 'G':
						*out_dinuc = dCG;
					break;
					case 'C':
						*out_dinuc = dCC;
					break;
					case 'T':
						*out_dinuc = dCT;
					break;
			}
			break;
			case 'T':
			switch(B)
			{
					case 'A':
						*out_dinuc = dTA;
					break;
					case 'G':
						*out_dinuc = dTG;
					break;
					case 'C':
						*out_dinuc = dTC;
					break;
					case 'T':
						*out_dinuc = dTT;
					break;
			}
			break;
			case '_':
			case 'N':

			break;
	}

	return NO_ERROR;
}

/**
 *
 * FILE OPERATIONS
 *
 */

/**
 * Print to file data from recalibration.
 */
ERROR_CODE
recal_fprint_info(const recal_info_t *data, const char *path)
{
	FILE *fp;
	int i,j;
	int n_quals = data->num_quals;
	int n_cycles = data->num_cycles;
	int n_dinuc = data->num_dinuc;

	fp = fopen(path, "w+");

	//Info msg
	printf("\n----------------\nPrinting data on \"%s\" file...\n----------------\n", path);

	//Print general info
	fprintf(fp, "==============================\nGENERAL: \nTotal miss: %lu\nTotal bases: %lu\nTotal delta: %.2f\n", data->total_miss, data->total_bases, data->total_delta);

	//Print quality infos
	fprintf(fp, "==============================\nQUAL VECTOR:\n");
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%3d %8.2f %10u %6.2f \n", i + MIN_QUALITY, data->qual_miss[i] /*- SMOOTH_CONSTANT_MISS*/, data->qual_bases[i] /*- SMOOTH_CONSTANT_BASES*/, data->qual_delta[i] /*- SMOOTH_CONSTANT_MISS*/);
	}

	//Print cycle infos
	fprintf(fp, "==============================\nQUAL-CYCLE MISS MATRIX (%d Quals - %d Cycles):\n", data->num_quals, data->num_cycles);
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%d \t", i + MIN_QUALITY);
		for(j = 0; j < n_cycles; j++)
		{
			fprintf(fp, "%.2f \t", data->qual_cycle_miss[i * n_cycles + j] /*- SMOOTH_CONSTANT_MISS*/);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "==============================\nQUAL-CYCLES BASES MATRIX (%d Quals - %d Cycles):\n", data->num_quals, data->num_cycles);
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%d \t", i + MIN_QUALITY);
		for(j = 0; j < n_cycles; j++)
		{
			fprintf(fp, "%u \t", data->qual_cycle_bases[i * n_cycles + j] /*- SMOOTH_CONSTANT_BASES*/);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "==============================\nQUAL-CYCLE DELTA MATRIX (%d Quals - %d Cycles):\n", data->num_quals, data->num_cycles);
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%d \t", i + MIN_QUALITY);
		for(j = 0; j < n_cycles; j++)
		{
			fprintf(fp, "%.2f \t", data->qual_cycle_delta[i * n_cycles + j]);
		}
		fprintf(fp, "\n");
	}

	//Print dinuc infos
	fprintf(fp, "==============================\nQUAL-DINUC MISS MATRIX (%d Quals - %d Dinuc):\n", data->num_quals, data->num_dinuc);
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%d \t", i + MIN_QUALITY);
		for(j = 0; j < n_dinuc; j++)
		{
			fprintf(fp, "%.2f \t", data->qual_dinuc_miss[i * n_dinuc + j] /*- SMOOTH_CONSTANT_MISS*/);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "==============================\nQUAL-DINUC BASES MATRIX (%d Quals - %d Dinuc):\n", data->num_quals, data->num_dinuc);
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%d \t", i + MIN_QUALITY);
		for(j = 0; j < n_dinuc; j++)
		{
			fprintf(fp, "%u \t", data->qual_dinuc_bases[i * n_dinuc + j] /*- SMOOTH_CONSTANT_BASES*/);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "==============================\nQUAL-DINUC DELTA MATRIX (%d Quals - %d Dinuc):\n", data->num_quals, data->num_dinuc);
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%d \t", i + MIN_QUALITY);
		for(j = 0; j < n_dinuc; j++)
		{
			fprintf(fp, "%.2f \t", data->qual_dinuc_delta[i * n_dinuc + j]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);

	return NO_ERROR;
}

/**
 * Save to file recalibration data.
 */
ERROR_CODE
recal_save_recal_info(const recal_info_t *data, const char *path)
{
	FILE *fp;

	printf("\n----------------\nSaving recalibration data to \"%s\"\n----------------\n", path);

	fp = fopen(path, "w+");

	fwrite(data, sizeof(qual_t), 1, fp);
	fwrite(data, sizeof(uint32_t), 3, fp);

	//Save total counters
	fwrite(&data->total_miss, sizeof(error_t), 1, fp);
	fwrite(&data->total_bases, sizeof(uint32_t), 1, fp);
	fwrite(&data->total_delta, sizeof(delta_t), 1, fp);

	//Save qual counters
	fwrite(data->qual_miss, sizeof(error_t), data->num_quals, fp);
	fwrite(data->qual_bases, sizeof(uint32_t), data->num_quals, fp);
	fwrite(data->qual_delta, sizeof(delta_t), data->num_quals, fp);

	//Save cycle counters
	fwrite(data->qual_cycle_miss, sizeof(error_t), data->num_quals * data->num_cycles, fp);
	fwrite(data->qual_cycle_bases, sizeof(uint32_t), data->num_quals * data->num_cycles, fp);
	fwrite(data->qual_cycle_delta, sizeof(delta_t), data->num_quals * data->num_cycles, fp);

	//Save dinuc counters
	fwrite(data->qual_dinuc_miss, sizeof(error_t), data->num_quals * data->num_dinuc, fp);
	fwrite(data->qual_dinuc_bases, sizeof(uint32_t), data->num_quals * data->num_dinuc, fp);
	fwrite(data->qual_dinuc_delta, sizeof(delta_t), data->num_quals * data->num_dinuc, fp);

	fclose(fp);

	return NO_ERROR;
}

/**
 * Load from file recalibration data.
 */
ERROR_CODE
recal_load_recal_info(const char *path, recal_info_t *data)
{
	FILE *fp;

	printf("\n----------------\nLoading recalibration data \"%s\"\n----------------\n", path);

	fp = fopen(path, "r");

	fread(data, sizeof(qual_t), 1, fp);
	fread(data, sizeof(uint32_t), 3, fp);

	//Read total counters
	fread(&data->total_miss, sizeof(error_t), 1, fp);
	fread(&data->total_bases, sizeof(uint32_t), 1, fp);
	fread(&data->total_delta, sizeof(delta_t), 1, fp);

	//Read qual counters
	fread(data->qual_miss, sizeof(error_t), data->num_quals, fp);
	fread(data->qual_bases, sizeof(uint32_t), data->num_quals, fp);
	fread(data->qual_delta, sizeof(delta_t), data->num_quals, fp);

	//Read cycle counters
	fread(data->qual_cycle_miss, sizeof(error_t), data->num_quals * data->num_cycles, fp);
	fread(data->qual_cycle_bases, sizeof(uint32_t), data->num_quals * data->num_cycles, fp);
	fread(data->qual_cycle_delta, sizeof(delta_t), data->num_quals * data->num_cycles, fp);

	//Read dinuc counters
	fread(data->qual_dinuc_miss, sizeof(error_t), data->num_quals * data->num_dinuc, fp);
	fread(data->qual_dinuc_bases, sizeof(uint32_t), data->num_quals * data->num_dinuc, fp);
	fread(data->qual_dinuc_delta, sizeof(delta_t), data->num_quals * data->num_dinuc, fp);

	fclose(fp);

	return NO_ERROR;
}

/**
 * PRIVATE FUNCTIONS
 */
