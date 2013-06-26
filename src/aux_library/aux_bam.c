#include "aux_bam.h"

void compare_bams_qual(const char* bamPath0, const char* bamPath1, int cycles)
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
}

bam_header_t* create_empty_bam_header(int num_chroms) {

	int i;

	//Memory allocation for struct
	bam_header_t *bam_header = (bam_header_t *) calloc(1, sizeof(bam_header_t));

	//Create a header with chroms targets number
	bam_header->n_targets = num_chroms;
	bam_header->target_name = (char **) calloc(num_chroms, sizeof(char *));
	bam_header->target_len = (uint32_t*) calloc(num_chroms, sizeof(uint32_t));

	for(i = 0; i < num_chroms; i++)
	{
		bam_header->target_name[i] = strdup("chr1");
		bam_header->target_len[i] = strlen("chr1")+1;
	}

	bam_header->text = strdup("@PG\tID:HPG-RECALIBRATOR\tVN:0.1\n");
	bam_header->l_text = strlen(bam_header->text);

	return bam_header;
}
