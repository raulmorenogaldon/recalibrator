#include "recal_structs.h"

/**
 *
 * DATA MANAGEMENT
 *
 */

/**
 * Allocate new recalibration data.
 */
recal_info_t *recal_new_info(int cycles)
{
	recal_info_t *data;
	int vector_size = MAX_QUALITY - MIN_QUALITY;

	data = malloc(sizeof(recal_info_t)); //Memory allocation

	//Struct initialization
	data->min_qual = MIN_QUALITY;
	data->num_quals = vector_size;
	data->num_cycles = cycles;
	data->num_dinuc = NUM_DINUC;

	//Total counters
	data->total_miss = 0;
	data->total_bases = 0;
	data->total_delta = 0;

	//Quality vectors
	data->qual_miss = new_vector(vector_size, SMOOTH_CONSTANT_MISS);
	data->qual_bases = new_vector(vector_size, SMOOTH_CONSTANT_BASES);
	data->qual_delta = new_vector_d(vector_size, 0.0);

	//Qual-Cycle matrix
	data->qual_cycle_miss = new_vector(vector_size * cycles, SMOOTH_CONSTANT_MISS);
	data->qual_cycle_bases = new_vector(vector_size * cycles, SMOOTH_CONSTANT_BASES);
	data->qual_cycle_delta = new_vector_d(vector_size * cycles, 0.0);

	//Qual-Dinuc matrix
	data->qual_dinuc_miss = new_vector(vector_size * cycles, SMOOTH_CONSTANT_MISS);
	data->qual_dinuc_bases = new_vector(vector_size * cycles, SMOOTH_CONSTANT_BASES);
	data->qual_dinuc_delta = new_vector_d(vector_size * cycles, 0.0);

	return data;
}

/**
 * Free recalibration data.
 */
void recal_destroy_info(recal_info_t *data)
{
	//Free quality vector
	free(data->qual_miss);
	free(data->qual_bases);
	free(data->qual_delta);

	//Free quality-cycle matrix
	free(data->qual_cycle_miss);
	free(data->qual_cycle_bases);
	free(data->qual_cycle_delta);

	//Free quality-dinuc matrix
	free(data->qual_dinuc_miss);
	free(data->qual_dinuc_bases);
	free(data->qual_dinuc_delta);

	//Free struct
	free(data);
}

/**
 * Add recalibration data from one base.
 */
void recal_add_base(recal_info_t *data, int qual, int cycle, int dinuc, int miss)
{
	int qual_index = qual - data->min_qual - 1;
	int qual_cycle_index = qual_index * data->num_cycles + cycle;
	int qual_dinuc_index = qual_index * data->num_dinuc + dinuc;

	if(qual - 1 < MIN_QUALITY_TO_STAT)
		return;

	//Error check
	if(qual_index >= MAX_QUALITY - MIN_QUALITY || qual_index < 0)
	{
		printf("add_base: ERROR, qual must be positive and minor than MAX_QUALITY - MIN_QUALITY ==> Qual = %d\n", qual);
		return;
	}
	if(cycle < 0 || cycle > data->num_cycles - 1)
	{
		printf("add_base: ERROR, cycle must be positive and minor than NUM_CYCLES\n");
		return;
	}
	if(dinuc > NUM_DINUC - 1)
	{
		printf("add_base: ERROR, dinuc must be minor than NUM_DINUC %d \n", dinuc);
		return;
	}

	//Increase total counters
	if(miss)
	data->total_miss++;
	data->total_bases++;

	//Increase quality vector
	if(miss)
	data->qual_miss[qual_index]++;
	data->qual_bases[qual_index]++;

	//Increase quality-cycle matrix
	if(miss)
	data->qual_cycle_miss[qual_cycle_index]++;
	data->qual_cycle_bases[qual_cycle_index]++;

	//Increase quality-dinuc matrix (if not dinuc != -1)
	if(dinuc >= 0)
	{
		if(miss)
		data->qual_dinuc_miss[qual_dinuc_index]++;
		data->qual_dinuc_bases[qual_dinuc_index]++;
	}
	else
	{
		printf("ERROR: unrecognized dinuc Q: %d, C: %d, D: %d, M: %d\n", qual, cycle, dinuc, miss);
		getchar();
	}
}

/**
 * Add recalibration data from vector of bases
 */
void recal_add_base_v(recal_info_t *data, char *seq, char *quals, int init_cycle, int end_cycle, char *dinuc, char *misses)
{
	int i;
	int qual_index;
	int qual_cycle_index;
	int qual_dinuc_index;

	//Iterates cycles
	for(i = init_cycle; i <= end_cycle; i++)
	{
		#ifdef NOT_COUNT_NUCLEOTIDE_N
		if(seq[i] != 'N' && quals[i] - 1 >= MIN_QUALITY_TO_STAT)
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
			if(!misses[i])
			{
				data->total_miss++;
				data->qual_miss[qual_index]++;
				data->qual_cycle_miss[qual_cycle_index]++;
				if(dinuc >= 0)
				{
					data->qual_dinuc_miss[qual_dinuc_index]++;
				}
			}
		}
	}
}

/**
 * Compute deltas from bases and misses.
 */
void recal_calc_deltas(recal_info_t* data)
{
	double global_empirical, phred, err0;
	double r_empirical;
	int matrix_index;
	int i, j;

	//Time measures
	#ifdef D_TIME_DEBUG
		time_init_slot(D_SLOT_CALC_DELTAS, clock(), TIME_GLOBAL_STATS);
	#endif

	printf("Processing deltas...\n");

	//Global delta
	global_empirical = (double)(data->total_miss + 1) / (double)(data->total_bases + 1);
	phred = 0.0;
	for(i = 0; i < data->num_quals; i++)
	{
		if(data->min_qual + i != 0)
		{
			//err0 = 1.0 / (10.0 * log10(data->min_qual + i + 1));
			//err0 = 1.0 / ((double)(data->min_qual + i));
			err0 = Pvalue((double)(data->min_qual + i));
			phred += err0 * (double)data->qual_bases[i];
		}
	}
	phred = phred / (double)data->total_bases;
	data->total_delta = Qvalue(global_empirical) - Qvalue(phred);

	//Delta R
	for(i = 0; i < data->num_quals; i++)
	{
		if(data->qual_bases[i] != 0)
		{
			data->qual_delta[i] = Qvalue(((double)(data->qual_miss[i]) / (double)(data->qual_bases[i])))
				- (double)(i + data->min_qual) - data->total_delta;
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
				data->qual_cycle_delta[matrix_index] = Qvalue((double)(data->qual_cycle_miss[matrix_index]) / (double)(data->qual_cycle_bases[matrix_index]))
					- (double)(i + data->min_qual) - (data->total_delta + data->qual_delta[i]);
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
				data->qual_dinuc_delta[matrix_index] = Qvalue((double)(data->qual_dinuc_miss[matrix_index]) / (double)(data->qual_dinuc_bases[matrix_index]))
					- (double)(i + data->min_qual) - (data->total_delta + data->qual_delta[i]);
			}
		}
	}

	printf("Deltas processed.\n");

	#ifdef D_TIME_DEBUG
		time_set_slot(D_SLOT_CALC_DELTAS, clock(), TIME_GLOBAL_STATS);
	#endif
}

/**
 *
 * ENUMERATION FUNCTIONS
 *
 */

/**
 * Return dinucleotide enumeration from two bases.
 */
enum DINUC recal_get_dinuc(char A, char B)
{
	switch(A)
	{
			case 'A':
			switch(B)
			{
					case 'A':
						return dAA;
					break;
					case 'G':
						return dAG;
					break;
					case 'C':
						return dAC;
					break;
					case 'T':
						return dAT;
					break;
			}
			break;
			case 'G':
			switch(B)
			{
					case 'A':
						return dGA;
					break;
					case 'G':
						return dGG;
					break;
					case 'C':
						return dGC;
					break;
					case 'T':
						return dGT;
					break;
			}
			break;
			case 'C':
			switch(B)
			{
					case 'A':
						return dCA;
					break;
					case 'G':
						return dCG;
					break;
					case 'C':
						return dCC;
					break;
					case 'T':
						return dCT;
					break;
			}
			break;
			case 'T':
			switch(B)
			{
					case 'A':
						return dTA;
					break;
					case 'G':
						return dTG;
					break;
					case 'C':
						return dTC;
					break;
					case 'T':
						return dTT;
					break;
			}
			break;
			case '_':
			case 'N':
				return d_X;
			break;
	}

	return -1;
}

/**
 *
 * FILE OPERATIONS
 *
 */

/**
 * Print to file data from recalibration.
 */
void recal_fprint_info(recal_info_t *data, const char *path)
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
		fprintf(fp, "%d \t%u \t%u \t%.2f \n", i + MIN_QUALITY, data->qual_miss[i] - SMOOTH_CONSTANT_MISS, data->qual_bases[i] - SMOOTH_CONSTANT_BASES, data->qual_delta[i] - SMOOTH_CONSTANT_MISS);
	}

	//Print cycle infos
	fprintf(fp, "==============================\nQUAL-CYCLE MISS MATRIX (%d Quals - %d Cycles):\n", data->num_quals, data->num_cycles);
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%d \t", i + MIN_QUALITY);
		for(j = 0; j < n_cycles; j++)
		{
			fprintf(fp, "%u \t", data->qual_cycle_miss[i * n_cycles + j] - SMOOTH_CONSTANT_MISS);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "==============================\nQUAL-CYCLES BASES MATRIX (%d Quals - %d Cycles):\n", data->num_quals, data->num_cycles);
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%d \t", i + MIN_QUALITY);
		for(j = 0; j < n_cycles; j++)
		{
			fprintf(fp, "%u \t", data->qual_cycle_bases[i * n_cycles + j] - SMOOTH_CONSTANT_BASES);
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
			fprintf(fp, "%u \t", data->qual_dinuc_miss[i * n_dinuc + j] - SMOOTH_CONSTANT_MISS);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "==============================\nQUAL-DINUC BASES MATRIX (%d Quals - %d Dinuc):\n", data->num_quals, data->num_dinuc);
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%d \t", i + MIN_QUALITY);
		for(j = 0; j < n_dinuc; j++)
		{
			fprintf(fp, "%u \t", data->qual_dinuc_bases[i * n_dinuc + j] - SMOOTH_CONSTANT_BASES);
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
}

/**
 * Save to file recalibration data.
 */
void recal_save_recal_info(recal_info_t *data, const char *path)
{
	FILE *fp;

	printf("\n----------------\nSaving recalibration data to \"%s\"\n----------------\n", path);

	fp = fopen(path, "w+");

	fwrite(data, sizeof(unsigned int), 4, fp);

	//Save total counters
	fwrite(&data->total_miss, sizeof(unsigned int), 1, fp);
	fwrite(&data->total_bases, sizeof(unsigned int), 1, fp);
	fwrite(&data->total_delta, sizeof(double), 1, fp);

	//Save qual counters
	fwrite(data->qual_miss, sizeof(unsigned int), data->num_quals, fp);
	fwrite(data->qual_bases, sizeof(unsigned int), data->num_quals, fp);
	fwrite(data->qual_delta, sizeof(double), data->num_quals, fp);

	//Save cycle counters
	fwrite(data->qual_cycle_miss, sizeof(unsigned int), data->num_quals * data->num_cycles, fp);
	fwrite(data->qual_cycle_bases, sizeof(unsigned int), data->num_quals * data->num_cycles, fp);
	fwrite(data->qual_cycle_delta, sizeof(double), data->num_quals * data->num_cycles, fp);

	//Save dinuc counters
	fwrite(data->qual_dinuc_miss, sizeof(unsigned int), data->num_quals * data->num_dinuc, fp);
	fwrite(data->qual_dinuc_bases, sizeof(unsigned int), data->num_quals * data->num_dinuc, fp);
	fwrite(data->qual_dinuc_delta, sizeof(double), data->num_quals * data->num_dinuc, fp);

	fclose(fp);
}

/**
 * Load from file recalibration data.
 */
void recal_load_recal_info(const char *path, recal_info_t *data)
{
	FILE *fp;

	printf("\n----------------\nLoading recalibration data \"%s\"\n----------------\n", path);

	fp = fopen(path, "r");

	fread(data, sizeof(unsigned int), 4, fp);

	//Read total counters
	fread(&data->total_miss, sizeof(unsigned int), 1, fp);
	fread(&data->total_bases, sizeof(unsigned int), 1, fp);
	fread(&data->total_delta, sizeof(double), 1, fp);

	//Read qual counters
	fread(data->qual_miss, sizeof(unsigned int), data->num_quals, fp);
	fread(data->qual_bases, sizeof(unsigned int), data->num_quals, fp);
	fread(data->qual_delta, sizeof(double), data->num_quals, fp);

	//Read cycle counters
	fread(data->qual_cycle_miss, sizeof(unsigned int), data->num_quals * data->num_cycles, fp);
	fread(data->qual_cycle_bases, sizeof(unsigned int), data->num_quals * data->num_cycles, fp);
	fread(data->qual_cycle_delta, sizeof(double), data->num_quals * data->num_cycles, fp);

	//Read dinuc counters
	fread(data->qual_dinuc_miss, sizeof(unsigned int), data->num_quals * data->num_dinuc, fp);
	fread(data->qual_dinuc_bases, sizeof(unsigned int), data->num_quals * data->num_dinuc, fp);
	fread(data->qual_dinuc_delta, sizeof(double), data->num_quals * data->num_dinuc, fp);

	fclose(fp);
}
