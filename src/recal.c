#include "recal.h"

//Struct creation functions
recal_info_t *recal_new_info(int cycles)
{
	recal_info_t *data;
	int vector_size = MAX_QUALITY - MIN_QUALITY;
	
	//Memory allocation
	data = malloc(sizeof(recal_info_t));
	
	//Struct initialization
	data->min_qual = MIN_QUALITY;
	data->num_quals = vector_size;
	data->num_cycles = cycles;
	data->num_dinuc = NUM_DINUC;
	
	data->total_miss = 0;	//Total counters
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

void initialize_vector(unsigned int *vector, int size, int value)
{
	int i;
	
	for(i = 0; i < size; i++)
	{
		vector[i] = value;
	}
}

unsigned int *new_vector(int size, int value)
{
	unsigned int *vector;
	int i;
	
	vector = (unsigned int *)malloc(size * sizeof(unsigned int));
	
	for(i = 0; i < size; i++)
		vector[i] = value;
		
	printf("Created new vector with %d positions, total size %lu bytes\n", size, size * sizeof(unsigned int));
	
	return vector;
}

double *new_vector_d(int size, double value)
{
	double *vector;
	int i;
	
	vector = (double *)malloc(size * sizeof(double));
	
	for(i = 0; i < size; i++)
		vector[i] = value;
		
	printf("Created new vector with %d positions, total size %lu bytes\n", size, size * sizeof(double));
	
	return vector;
}

//Manage functions
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
		printf("ERROR: unrecognised dinuc Q: %d, C: %d, D: %d, M: %d\n", qual, cycle, dinuc, miss);
		getchar();
	}
}

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

//Preprocess fasta reference
void recal_preprocess_fasta(char* input_fasta, char* output_file)
{
	generate_codes(output_file, input_fasta);
}

//Obtain data from BAM
void recal_get_data_from_file(char *bam_path, char *ref_name, char *ref_path, recal_info_t *out_info)
{
	genome_t* ref;
	bam_file_t *bam_f;
	
	//Open bam
	printf("Opening BAM from \"%s\" ...\n", bam_path);
	bam_f = bam_fopen(bam_path);
	printf("BAM opened!...\n");
	
	//Open reference genome
	printf("Opening reference genome from \"%s%s\" ...\n", ref_path, ref_name);	
	ref = genome_new(ref_name, ref_path);
	printf("Reference opened!...\n");
	
	//Fill data
	#ifdef D_TIME_DEBUG
		time_init_slot(D_SLOT_GET_DATA_BAM, clock(), time_global_stats);
	#endif
	recal_get_data_from_bam(bam_f, ref, out_info);
	#ifdef D_TIME_DEBUG
		time_set_slot(D_SLOT_GET_DATA_BAM, clock(), time_global_stats);
	#endif
	
	//Memory free	
	printf("\nClosing BAM file...\n");
	bam_fclose(bam_f);
	genome_free(ref);
	printf("BAM closed.\n");
}

void recal_get_data_from_bam(bam_file_t *bam, genome_t* ref, recal_info_t* output_data)
{
	bam_batch_t* batch;
	int count = 0;
	
	#ifdef USE_BATCH_POOL
	bam_pool_t* pool;
	#endif
	
	int zero = 0;
	
	unmapped = 0;
	
	//DEBUG
	#ifdef D_SAMPLES_OUTPUT
	FILE *fp;
	out_path = "cadenas.out";
	fp = fopen(out_path, "w");
	fclose(fp);
	#endif
	
	printf("\n----------------\nProcessing \"%s\" file...\n----------------", bam->filename);
	
	//Allocate memory for batchs
	batch = bam_batch_new(MAX_BATCH_SIZE, MULTIPLE_CHROM_BATCH);
	
	//Read batch
	#ifdef D_TIME_DEBUG
		time_init_slot(D_SLOT_READ_BATCH, clock(), time_global_stats);
	#endif
	bam_fread_max_size(batch, MAX_BATCH_SIZE, 1, bam);
	#ifdef D_TIME_DEBUG
		time_set_slot(D_SLOT_READ_BATCH, clock(), time_global_stats);
	#endif
		
	printf("\nNum alignments in batchs: %d\n----------------\n", batch->num_alignments);
	
	#ifdef USE_BATCH_POOL
	printf("Using batch pool!\n");
	pool = pool_new(3, MAX_BATCH_SIZE);
	#endif
	
	while(batch->num_alignments != 0 
		#ifdef D_MAX_READS 
			&& count < D_MAX_READS 
		#endif
		)
	{
		#ifdef USE_BATCH_POOL
		while(!pool_has_free_batchs(pool));
		{
			#ifdef D_TIME_DEBUG
				time_init_slot(D_SLOT_MEMCOPY_BATCH, clock(), time_global_stats);
			#endif
			pool_add_batch(batch, pool);
			bam_batch_free(batch, 1);
			
			//Get batch in pool
			pool_next_batch(pool);
			batch = pool_get_batch(pool);
			#ifdef D_TIME_DEBUG
				time_set_slot(D_SLOT_MEMCOPY_BATCH, clock(), time_global_stats);
			#endif
		}
		#endif
		
		//Process batch
		#ifdef D_TIME_DEBUG
			time_init_slot(D_SLOT_GET_DATA_BATCH, clock(), time_global_stats);
		#endif
		recal_get_data_from_bam_batch(batch, ref, output_data);
		#ifdef D_TIME_DEBUG
			time_set_slot(D_SLOT_GET_DATA_BATCH, clock(), time_global_stats);
		#endif
	
		//Update read counter
		count += batch->num_alignments;
		
		//Show total progress
		printf("\n Total alignments readed: %d", count);	
		
		//Free memory and take a new batch
		#ifndef USE_BATCH_POOL
		bam_batch_free(batch, 1);
		#endif
		batch = bam_batch_new(MAX_BATCH_SIZE, MULTIPLE_CHROM_BATCH);
		
		//Read batch
		#ifdef D_TIME_DEBUG
			time_init_slot(D_SLOT_READ_BATCH, clock(), time_global_stats);
		#endif
		bam_fread_max_size(batch, MAX_BATCH_SIZE, 1, bam);	
		#ifdef D_TIME_DEBUG
			time_set_slot(D_SLOT_READ_BATCH, clock(), time_global_stats);
		#endif
	} 
	
	printf("\n----------------\n%d alignments readed.", count);
	printf("\n%d alignments processed.", count - unmapped);
	printf("\n%d alignments unmapped.", unmapped);
		
	#ifdef USE_BATCH_POOL
		pool_free(pool);
	#endif
	
	//Free batch
	bam_batch_free(batch, 1);
}

void recal_get_data_from_bam_batch(bam_batch_t* batch, genome_t* ref, recal_info_t* output_data)
{
	int i;
	
	//Process all alignments of the batchs
	for(i = 0; i < batch->num_alignments; i++)
	{
		//Process every alignment
		#ifdef NOT_MAPPING_QUAL_ZERO
		if(batch->alignments_p[i]->core.qual != 0)
		#endif
		{
		#ifdef D_TIME_DEBUG
			time_init_slot(D_SLOT_GET_DATA_ALIG, clock(), time_global_stats);
		#endif
		recal_get_data_from_bam_alignment(batch->alignments_p[i], ref, output_data); 
		#ifdef D_TIME_DEBUG
			time_set_slot(D_SLOT_GET_DATA_ALIG, clock(), time_global_stats);
		#endif
		}
	}
}

inline void recal_get_data_from_bam_alignment(bam1_t* alig, genome_t* ref, recal_info_t* output_data)
{
	//alignment_t* aux_alig;
	FILE *fp;
	char *quals;
	char *ref_seq;
	char *bam_seq, *bam_seq_aux;
	char aux_comp[16];
	char *comp_res, *dinucs;
	int i, miss, j;
	#ifdef __SSE2__
	__m128i v_ref, v_seq, v_comp;
	#endif
	enum DINUC dinuc;
	unsigned long int init_pos, end_pos;
	uint32_t flag; 
	
	//Convert bam1_t to alignment_t
	//aux_alig = alignment_new_by_bam(alig, 1);
	
	//Unmapped readings	
	if(alig->core.tid < 0)
	{
		unmapped++;
		return;
	}

	//Bam seq fields
	bam_seq = convert_to_sequence_string(bam1_seq(alig), alig->core.l_qseq);

	//Qual fields
	quals = (char*) calloc(alig->core.l_qseq + 1, sizeof(char));
	convert_to_quality_string_length(quals, bam1_qual(alig), alig->core.l_qseq, 1);

	//Sequence fields
	ref_seq = malloc(alig->core.l_qseq * sizeof(char));
	flag = (uint32_t) alig->core.flag;
	init_pos = alig->core.pos + 1;
	end_pos = init_pos + (alig->core.l_qseq - 1);
	
	//Obtain reference for this 100 nucleotides
	genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, alig->core.tid, &init_pos, &end_pos, ref);
	
	//DEBUG
	#ifdef D_SAMPLES_OUTPUT
	if(rand() % D_PROB_SAMPLE == 1)
	{
		fp = fopen(out_path, "a+");
		fprintf(fp, "--------------------------\n");
		fprintf(fp, "Chrom:%d\nPos ref:%lu\nPos bam:%d\nL_qseq:%d\nSEQ:%s\nBAM:%s\n", alig->core.tid, init_pos, alig->core.pos, alig->core.l_qseq, ref_seq, bam_seq);
		fprintf(fp, "Mis:");
		for(i = 0; i < alig->core.l_qseq; i++)
		{
			fprintf(fp, "%d",  ref_seq[i] != bam_seq[i]);
		}
		fprintf(fp,"\n");
		fclose(fp);
	}
	#endif
	
	//Iterates nucleotides in this read
	dinuc = 0;
	comp_res = malloc(alig->core.l_qseq * sizeof(char));
	dinucs = malloc(alig->core.l_qseq * sizeof(int));
	for(i = 0; i < alig->core.l_qseq; i++)
	{	
		#ifdef __SSE2__ /*SSE Block*/
		if( (i + 16) < alig->core.l_qseq)
		{
			//Use SSE
			//_mm_prefetch(&ref_seq[i + 16], _MM_HINT_T0);
			//_mm_prefetch(&bam_seq[i + 16], _MM_HINT_T0);

			//Pack sequences
			v_ref = _mm_load_si128(&ref_seq[i]);
			v_seq = _mm_load_si128(&bam_seq[i]);
			
			//Compare sequences
			v_comp = _mm_cmpeq_epi8(v_ref, v_seq);
			
			//Store comparation values
			_mm_store_si128(&comp_res[i], v_comp);
			
			/*for(j = i; j < i + 16; j++)
			{
				#ifdef NOT_COUNT_NUCLEOTIDE_N
				if(bam_seq[j] != 'N')
				#endif
				{
					//Only take dinuc if we arent in first cycle
					if(j > 0) 
					{
						dinuc = recal_get_dinuc(bam_seq[j-1], bam_seq[j]);
					} 
					else
					{
						dinuc = recal_get_dinuc('_', '_');
					}
					
					//if this happen, then error will be thrown
					if(dinuc == -1)
					{
						dinuc = d_X;
						//printf("\nDinuc:%c%c", bam_seq[i-1],bam_seq[i]);
					}
					
					//Add data
					recal_add_base(output_data, quals[j], j, dinuc, v_comp[j-i]);
				}
			}*/
			i += 15;
		}
		else
		#endif /*SSE Block*/
		{
			/*#ifdef NOT_COUNT_NUCLEOTIDE_N
			if(bam_seq[i] != 'N')
			#endif
			{
				//Normal
				//Only take dinuc if we arent in first cycle
				if(i > 0) 
				{
					dinuc = recal_get_dinuc(bam_seq[i-1], bam_seq[i]);
				} 
				else
				{
					dinuc = recal_get_dinuc('_', '_');
				}
				
				//if this happen, then error will be thrown
				if(dinuc == -1)
				{
					dinuc = d_X;
					//printf("\nDinuc:%c%c", bam_seq[i-1],bam_seq[i]);
				}
					
				//Add data
				recal_add_base(output_data, quals[i], i, dinuc, ref_seq[i] != bam_seq[i]);
			}*/
			if(ref_seq[i] != bam_seq[i])
			{
				comp_res[i] = 0;	/*Diff*/
			}
			else
			{
				comp_res[i] = -1;	/*Equals*/
			}
		}		
	}
	
	//Dinucs
	for(i = 0; i < alig->core.l_qseq; i++)
	{
		if(i > 0) 
		{
			dinucs[i] = recal_get_dinuc(bam_seq[i-1], bam_seq[i]);
			if(dinucs[i] == -1)
			{
				dinucs[i] = d_X;
			}
		} 
		else
		{
			dinucs[i] = d_X;
		}
	}
	
	//Add data
	recal_add_base_v(output_data, bam_seq, quals, 0, alig->core.l_qseq - 1, dinucs, comp_res);

	free(comp_res);
	free(dinucs);
	free(quals);
	free(ref_seq);
	free(bam_seq);
	
	//alignment_free(aux_alig);
}

void recal_calc_deltas(recal_info_t* data)
{
	double global_empirical, phred, err0;
	double r_empirical;
	int matrix_index;
	int i, j;
	
	//Time measures
	#ifdef D_TIME_DEBUG
		time_init_slot(D_SLOT_CALC_DELTAS, clock(), time_global_stats);
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
		time_set_slot(D_SLOT_CALC_DELTAS, clock(), time_global_stats);
	#endif
}

//Recalibrate bam function
void recal_recalibrate_bam_file(char *orig_bam_path, recal_info_t *bam_info, char *recal_bam_path)
{
	bam_file_t *orig_bam_f, *recal_bam_f;
	bam_header_t *recal_bam_header;
	
	//Open bam
	printf("Opening BAM from \"%s\" to being recalibrated ...\n", orig_bam_path);
	orig_bam_f = bam_fopen(orig_bam_path);
	printf("BAM opened!...\n");
	
	//Create new bam
	printf("Creating new bam file in \"%s\"...\n", recal_bam_path);
	recal_bam_header = create_empty_bam_header(orig_bam_f->bam_header_p->n_targets);
	recal_bam_f = bam_fopen_mode(recal_bam_path, recal_bam_header, "w");
	bam_fwrite_header(recal_bam_header, recal_bam_f);
	printf("New BAM initialized!...\n");
	
	//Recalibrate bams
	#ifdef D_TIME_DEBUG
		time_init_slot(D_SLOT_RECALIBRATE, clock(), time_global_stats);
	#endif
	recal_recalibrate_bam(orig_bam_f, bam_info, recal_bam_f);
	#ifdef D_TIME_DEBUG
		time_set_slot(D_SLOT_RECALIBRATE, clock(), time_global_stats);
	#endif
	
	//Memory free	
	printf("Closing \"%s\" BAM file...\n", orig_bam_path);
	bam_fclose(orig_bam_f);
	printf("Closing \"%s\" BAM file...\n", recal_bam_path);
	bam_fclose(recal_bam_f);
	printf("BAMs closed.\n");
	printf("Recalibration DONE.\n");
}

void recal_recalibrate_bam(bam_file_t *orig_bam_f, recal_info_t *bam_info, bam_file_t *recal_bam_f)
{
	bam_batch_t* batch;
	int count = 0;
	
	#ifdef USE_BATCH_POOL
	bam_pool_t* pool;
	#endif
	
	//Allocate memory for batchs
	batch = bam_batch_new(MAX_BATCH_SIZE, MULTIPLE_CHROM_BATCH);
	
	//Read batch
	bam_fread_max_size(batch, MAX_BATCH_SIZE, 1, orig_bam_f);
	
	#ifdef USE_BATCH_POOL
	printf("Using batch pool!\n");
	pool = pool_new(3, MAX_BATCH_SIZE);
	#endif
	
	#ifdef D_RECAL_INTER_RESULTS
		printf("=======================\nQual intermediate results\n");
	#endif
	
	while(batch->num_alignments != 0
		#ifdef D_MAX_READS_W 
			&& count < D_MAX_READS_W 
		#endif
		)
	{
		#ifdef USE_BATCH_POOL
		if(pool_has_free_batchs(pool))
		{
			pool_add_batch(batch, pool);
			bam_batch_free(batch, 1);
			
			//Get batch in pool
			pool_next_batch(pool);
			batch = pool_get_batch(pool);
		}
		#endif
		
		//Recalibrate batch
		recal_recalibrate_batch(batch, bam_info, recal_bam_f);

		//Update read counter
		count += batch->num_alignments;
		
		//Show total progress
		printf("\n\tTotal alignments recalibrated: %d\n", count);	
				
		//Free memory and take a new batch
		#ifndef USE_BATCH_POOL
		bam_batch_free(batch, 1);
		#endif
		batch = bam_batch_new(MAX_BATCH_SIZE, MULTIPLE_CHROM_BATCH);
		
		//Read batch
		bam_fread_max_size(batch, MAX_BATCH_SIZE, 1, orig_bam_f);	
		
	} 
	#ifdef D_RECAL_INTER_RESULTS
		printf("\n=======================\n");
	#endif
	
	#ifdef USE_BATCH_POOL
		pool_free(pool);
	#endif
	
	printf("\n---------------------\n", count);
}

void recal_recalibrate_batch(bam_batch_t* batch, recal_info_t *bam_info, bam_file_t *recal_bam_f)
{
	int i;
	
	//Process all alignments of the batchs
	for(i = 0; i < batch->num_alignments; i++)
	{
		#ifdef D_TIME_DEBUG
		time_init_slot(D_SLOT_RECAL_ALIG, clock(), time_global_stats);
		#endif
		//Process every alignment
		recal_recalibrate_alignment(batch->alignments_p[i], bam_info, recal_bam_f); 
		#ifdef D_TIME_DEBUG
		time_set_slot(D_SLOT_RECAL_ALIG, clock(), time_global_stats);
		#endif
	}
}

void recal_recalibrate_alignment(bam1_t* alig, recal_info_t *bam_info, bam_file_t *recal_bam_f)
{
	int qual_index;
	int matrix_index;
	int i;
	enum DINUC dinuc;
	double delta_r, delta_rc, delta_rd;
	char *quals;
	char *bam_seq;
	
	alignment_t* aux_alig;
	bam1_t *aux_alig1;
	
	//Convert bam1_t to alignment_t
	aux_alig = alignment_new_by_bam(alig, 1);
	
	//Bam seq fields
	bam_seq = convert_to_sequence_string(bam1_seq(alig), alig->core.l_qseq);
	
	//Qual fields
	quals = (char*) calloc(alig->core.l_qseq + 1, sizeof(char));
	convert_to_quality_string_length(quals, bam1_qual(alig), alig->core.l_qseq, 1);

	//Iterates nucleotides in this read
	dinuc = 0;
	for(i = 0; i < alig->core.l_qseq; i++)
	{			
		//Compare only if the nucleotide is not "N"
		#ifdef NOT_COUNT_NUCLEOTIDE_N
		if(bam_seq[i] != 'N')
		#endif
		{
			//Recalibrate quality
			qual_index = quals[i] - bam_info->min_qual;
			delta_r = bam_info->qual_delta[qual_index];
			
			matrix_index = qual_index * bam_info->num_cycles + i;
			delta_rc = bam_info->qual_cycle_delta[matrix_index];
			
			//dont take prev dinuc in first cycle (delta = 0)
			if(i > 0) 
			{
				dinuc = recal_get_dinuc(bam_seq[i-1], bam_seq[i]);
				matrix_index = qual_index * bam_info->num_dinuc + i;
				delta_rd = bam_info->qual_dinuc_delta[matrix_index];
			} 
			else
			{
				delta_rd = 0.0;
			}
			
			#ifdef D_RECAL_INTER_RESULTS
				int aux_qual = quals[i];
			#endif
			
			//Recalibration formula
			quals[i] = quals[i] + bam_info->total_delta + delta_r + delta_rc + delta_rd;
			
			#ifdef D_RECAL_INTER_RESULTS
			if(rand() % D_RECAL_INTER_PROB == 1)
			{
				printf("Q: %d ",aux_qual);
				printf("---- %d\n",quals[i]);
				printf("Cy: %d, Qr + DQg:%f + DQr:%f + DQc:%f + DQd:%f\n", i, bam_info->total_delta, delta_r, delta_rc, delta_rd);
			}
			#endif
		}
	}
	
	//Set qualities in alignment
	free(aux_alig->quality);
	aux_alig->quality = quals;
	aux_alig1 = convert_to_bam(aux_alig, 1);
	
	//Write in bam
	bam_fwrite(aux_alig1, recal_bam_f);
	
	//Memory free
	alignment_free(aux_alig);
	bam_destroy1(aux_alig1);
	free(bam_seq);
}

//Get dinuc
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

//Q and P
inline double Qvalue(double P)
{
	#ifdef P_SOLEXA
		return Qsolexa(P);
	#else
		return Qsanger(P);
	#endif
}

double Pvalue(double Q)
{
	#ifdef P_SOLEXA
		return Psolexa(Q);
	#else
		return Psanger(Q);
	#endif
}

//Free functions
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

//File print info
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

//File recal_info save
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
