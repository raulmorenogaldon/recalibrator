#include "aux_bam.h"

ERROR_CODE
compare_bams_qual(const char* bamPath0, const char* bamPath1, const int cycles)
{
	bam_file_t* bamFile0;
	bam_file_t* bamFile1;
	bam_batch_t* bamBatch0;
	bam_batch_t* bamBatch1;
	bam1_t* bamAlig;
	alignment_t* aligAlig0;
	alignment_t* aligAlig1;
	int diff, i;

	printf("Opening BAM 1 form \"%s\" ...\n", bamPath0);
	printf("Opening BAM 2 form \"%s\" ...\n", bamPath1);
	bamFile0 = bam_fopen(bamPath0);
	bamFile1 = bam_fopen(bamPath1);
	printf("BAM opened!...\n");


	printf("\n\n---------------------------------------------------------\n");

	bamBatch0 = bam_batch_new(1, MULTIPLE_CHROM_BATCH);
	bamBatch1 = bam_batch_new(1, MULTIPLE_CHROM_BATCH);
	bam_fread_max_size(bamBatch0, 1, 1, bamFile0);
	bam_fread_max_size(bamBatch1, 1, 1, bamFile1);

	//Obtain first alignment from first bam
	bamAlig = bamBatch0->alignments_p[0];
	aligAlig0 = alignment_new_by_bam(bamAlig, 1);
	alignment_print(aligAlig0);

	//Obtain first alignment from second bam
	bamAlig = bamBatch1->alignments_p[0];
	aligAlig1 = alignment_new_by_bam(bamAlig, 1);
	alignment_print(aligAlig1);

	//Obtain quality diffs
	printf("Diffs: \nNuc\tQ1\tQ2\n");
	diff=0;
	for(i=0; i < 76; i++)
	{
		printf("%c \t%d ", aligAlig0->sequence[i], aligAlig0->quality[i]);
		if(aligAlig0->quality[i] == aligAlig1->quality[i])
		{
			printf("====\t%d\n", aligAlig1->quality[i]);
		}
		else
		{
			printf("\t%d\n", aligAlig1->quality[i]);
		}


		diff += abs(aligAlig1->quality[i] - aligAlig0->quality[i]);
	}
	printf("Total diff: %d\n", diff);

	printf("\n---------------------------------------------------------\n");
	printf("Closing BAMs...\n");
	bam_fclose(bamFile0);
	bam_fclose(bamFile1);
	bam_batch_free(bamBatch0, 1);
	bam_batch_free(bamBatch1, 1);

	printf("BAM closed.\n");

	return NO_ERROR;
}

ERROR_CODE
init_empty_bam_header(const unsigned int num_chroms, bam_header_t *header)
{
	int i;

	if(num_chroms == 0)
		return INVALID_INPUT_PARAMS_0;

	if(!header)
	{
		return INVALID_INPUT_PARAMS_NULL;
	}

	//Create a header with chroms targets number
	header->n_targets = num_chroms;
	header->target_name = (char **) calloc(num_chroms, sizeof(char *));
	header->target_len = (uint32_t*) calloc(num_chroms, sizeof(uint32_t));

	for(i = 0; i < num_chroms; i++)
	{
		header->target_name[i] = strdup("chr1");
		header->target_len[i] = strlen("chr1")+1;
	}

	header->text = strdup("@PG\tID:HPG-RECALIBRATOR\tVN:0.1\n");
	header->l_text = strlen(header->text);

	return NO_ERROR;
}

char *
new_sequence_from_bam(bam1_t *bam1)
{
	char *bam_seq = bam1_seq(bam1);
	int seq_len = bam1->core.l_qseq;

	char *seq = (char *) _mm_malloc(seq_len * sizeof(char), MEM_ALIG_SIZE);

	// nucleotide content
	for (int i = 0; i < seq_len; i++) {
		switch (bam1_seqi(bam_seq, i))
		{
		case 1:
			seq[i] = 'A';
			break;
		case 2:
			seq[i] = 'C';
			break;
		case 4:
			seq[i] = 'G';
			break;
		case 8:
			seq[i] = 'T';
			break;
		case 15:
			seq[i] = 'N';
			//printf("N");
			break;
		default:
			seq[i] = 'N';
			break;
		}
	}

	return seq;
}

char *
new_quality_from_bam(bam1_t *bam1, int base_quality)
{
	char *bam_qual = bam1_qual(bam1);
	int qual_len = bam1->core.l_qseq;

	char *qual = (char *) malloc(qual_len * sizeof(char));
	for (int i = 0; i < qual_len; i++) {
		qual[i] = base_quality + bam_qual[i];
	}

	return qual;
}

ERROR_CODE
decompose_cigar(char *cigar, uint8_t cigar_l, char *n_elem, char *type, uint8_t *types_l, uint8_t max_types_length)
{
	uint8_t u_elems;
	char c_type;
	char cigar_elem;
	int pos;
	int i;

	if(cigar == NULL
			|| cigar_l <= 0
			|| n_elem == NULL
			|| type == NULL
			|| types_l == NULL)
	{
		return INVALID_INPUT_PARAMS_NULL;
	}

	pos = 0;
	n_elem[pos] = 0;
	for(i = 0; i < cigar_l; i++)
	{
		cigar_elem = cigar[i];

		//Get number of elem
		if(cigar_elem <= '9' && cigar_elem >= '0')
		{
			//Is number
			n_elem[pos] *= 10;
			n_elem[pos] += cigar_elem - '0';
		}
		else
		{
			//Is type
			//switch(cigar_elem)
			//{
			//case 'I':
			//case 'D':
			//default:	//Missmatch, etc...
				if(n_elem[pos] > 0)	//Avoid void types like '0D', '0M' ..
				{
					type[pos] = cigar_elem;
					pos++;

					if(pos >= max_types_length)
					{
						//Reached max number of cigar components
						goto decompose_cigar_end;
					}

					n_elem[pos] = 0;
				}
				//break;
			//}
		}
	}

	//Set type length
	decompose_cigar_end:
	if(types_l != NULL)
		*types_l = pos;

	return NO_ERROR;
}

ERROR_CODE
supress_indels(char *seq, uint8_t seq_l, char *cigar_elem, char *cigar_type, uint8_t cigar_type_l, char *seq_res, uint8_t *seq_res_l)
{
	int i, j, seq_i, res_i;
	char count;

	//Check null parameters
	if(seq == NULL
		|| cigar_elem == NULL
		|| cigar_type == NULL
		|| cigar_type_l == NULL
		|| seq_res == NULL
		|| seq_res_l == NULL)
	{
		return INVALID_INPUT_PARAMS_NULL;
	}

	//Iterate cigar
	for(i = 0, res_i = 0, seq_i = 0; i < cigar_type_l; i++)
	{
		switch(cigar_type[i])
		{
		case 'I':	//Insertion
			seq_i += cigar_elem[i];
			break;
		case 'D':	//Deletion
			count = cigar_elem[i];
			for(j = 0; j < count; j++)
			{
				seq_res[res_i] = 'X';
				res_i++;
			}
			break;
		default:	//Missmatch, etc...
			count = cigar_elem[i];
			memcpy(&seq_res[res_i], &seq[seq_i], count);
			res_i += count;
			seq_i += count;
			break;
		}
	}

	//Ser result length
	*seq_res_l = res_i;

	//Set string null character in las position
	seq_res[res_i] = '\0';

	return NO_ERROR;
}

ERROR_CODE
supress_indels_from_32_cigar(char *seq, char *qual, uint8_t seq_l, uint32_t *cigar, uint8_t cigar_l, char *seq_res, char *qual_res, uint8_t *seq_res_l)
{
	int i, j, seq_i, res_i;
	char count;
	uint32_t elem;
	uint8_t type;

	//Check null parameters
	if(seq == NULL
		|| cigar == NULL
		|| seq_res == NULL
		|| seq_res_l == NULL)
	{
		return INVALID_INPUT_PARAMS_NULL;
	}

	for(i = 0, seq_i = 0, res_i = 0; i < cigar_l; i++)
	{
		elem = cigar[i];
		count = elem >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		type = elem & BAM_CIGAR_MASK;	//Get type from cigar

		switch(type)
		{
		case BAM_CINS:	//Insertion
			seq_i += count;
			break;
		case BAM_CDEL:	//Deletion
			for(j = 0; j < count; j++)
			{
				seq_res[res_i] = 'X';
				if(qual && qual_res) qual_res[res_i] = '!';
				res_i++;
			}
			break;
		default:	//Missmatch, etc...
			memcpy(&seq_res[res_i], &seq[seq_i], count);
			if(qual && qual_res) memcpy(&qual_res[res_i], &qual[seq_i], count);
			res_i += count;
			seq_i += count;
			break;
		}
	}

	//Ser result length
	*seq_res_l = res_i;

	//Set string null character in last position
	seq_res[res_i] = '\0';
	qual_res[res_i] = '\0';

	return NO_ERROR;
}
