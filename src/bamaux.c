#include "bamaux.h"

//Struct creation
bam_info_t *bam_new_info()
{
	bam_info_t *info;
	
	info = malloc(sizeof(bam_info_t));
	
	//Field init to 0
	info->min_qual = 0;
	info->max_qual = 0;
	info->num_cycles = 0;
	
	return info;
}

//BAM analysis
void bam_get_info(bam_file_t *bam, bam_info_t *info)
{
	bam_batch_t* batch;
	int i, count;
	
	printf("\n----------------\nGathering info from \"%s\" file...\n----------------", bam->filename);
	
	//Allocate memory for batchs
	batch = bam_batch_new(INFO_BATCH_SIZE, MULTIPLE_CHROM_BATCH);
	
	//Read batch
	bam_fread_max_size(batch, INFO_BATCH_SIZE, 1, bam);
	
	count = 0;
	while(batch->num_alignments != 0)
	{
		//Process all alignments of the batchs
		for(i = 0; i < batch->num_alignments; i++)
		{
			//Process every alignment
			bam_get_info_from_alignment(batch->alignments_p[i], info); 
		}
		
		//Show total progress
		count += batch->num_alignments;
		printf("\nTotal alignments readed: %d", count);	
		
		//Free memory and take a new batch
		bam_batch_free(batch, 1);
		batch = bam_batch_new(INFO_BATCH_SIZE, MULTIPLE_CHROM_BATCH);
		
		//Read batch
		bam_fread_max_size(batch, INFO_BATCH_SIZE, 1, bam);	
	} 
	
	printf("\n----------------\n%d alignments readed.", count);
	
	//Free batch
	bam_batch_free(batch, 1);
}

void bam_get_info_from_alignment(bam1_t* alig, bam_info_t* info)
{
	char *quals;
	int i;
	
	//Get num cycles
	if(info->num_cycles < alig->core.l_qseq)
	{
		info->num_cycles = alig->core.l_qseq;
	}
	
	//Qual fields
	quals = (char*) calloc(alig->core.l_qseq + 1, sizeof(char));
	convert_to_quality_string_length(quals, bam1_qual(alig), alig->core.l_qseq, 1);

	//Update quals limits
	for(i = 0; i < alig->core.l_qseq; i++)
	{
		if(quals[i] < info->min_qual)
		{
			info->min_qual = quals[i];
		}
		else 
		{
			if(quals[i] > info->max_qual)
			{
				info->max_qual = quals[i];
			}
		}
	}
	
	free(quals);
}

//Struct destroy
void bam_destroy_info(bam_info_t *info)
{
	free(info);
}

//Create an empty bam header
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

//Math operations
inline double Qsolexa(double p)
{
	return -10.0 * log10(p/(1-p));
}

inline double Psolexa(double Q)
{
	return pow(10.0, -Q/10.0) / (pow(10.0, -Q/10.0) + 1);
}

double Qsanger(double p)
{
	return -10.0 * log10(p);
}

inline double Psanger(double Q)
{
	return pow(10.0, -Q/10.0);
}

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

//Optimization routines
#define cpuid(func,ax,bx,cx,dx)\
	__asm__ __volatile__ ("cpuid":\
	"=a" (ax), "=b" (bx), "=c" (cx), "=d" (dx) : "a" (func));

void printf_proc_features()
{
	/*int i;
	char vector[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
	char vector2[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
	__m128i v, v2;
	
	for(i = 0; i < 1000000000; i++)
	{
		
	v = _mm_load_si128(vector);
	v2 = _mm_load_si128(vector2);
	
	v = _mm_add_epi8(v, v2);
	
	}
	
	_mm_stream_si128(vector, v);
	
	printf("Vector: ");
	for(i = 0; i < 16; i++)
	{
		printf("%d ",vector[i]);
	}
	printf("\n");*/
	
	int x64     = 0;
	int MMX     = 0;
	int SSE     = 0;
	int SSE2    = 0;
	int SSE3    = 0;
	int SSSE3   = 0;
	int SSE41   = 0;
	int SSE42   = 0;
	int SSE4a   = 0;
	int AVX     = 0;
	int XOP     = 0;
	int FMA3    = 0;
	int FMA4    = 0;
	int a,b,c,d;
	
	cpuid(0,a,b,c,d);
	int nIds = a;
	
	cpuid(0x80000000,a,b,c,d);
	int nExIds = a;
	
	if (nIds >= 1)
	{
		cpuid(0x00000001,a,b,c,d);
		MMX   = (d & ((int)1 << 23)) != 0;
		SSE   = (d & ((int)1 << 25)) != 0;
		SSE2  = (d & ((int)1 << 26)) != 0;
		SSE3  = (c & ((int)1 <<  0)) != 0;

		SSSE3 = (c & ((int)1 <<  9)) != 0;
		SSE41 = (c & ((int)1 << 19)) != 0;
		SSE42 = (c & ((int)1 << 20)) != 0;

		AVX   = (c & ((int)1 << 28)) != 0;
		FMA3  = (c & ((int)1 << 12)) != 0;
    }
    
    if (nExIds >= 0x80000001){
		cpuid(0x80000001,a,b,c,d);
    
		x64   = (d & ((int)1 << 29)) != 0;
		SSE4a = (c & ((int)1 <<  6)) != 0;
		FMA4  = (c & ((int)1 << 16)) != 0;
		XOP   = (c & ((int)1 << 11)) != 0;
	}
	
	printf("-----------------\nProcessor features:\n");
	
	if(x64)
		printf("64 bits\n");
	
	printf("Supported: ");
	
	if(MMX)
		printf("/MMX");
		
	if(SSE)
		printf("/SSE");
		
	if(SSE2)
		printf("/SSE2");
	
	if(SSE3)
		printf("/SSE3");
		
	printf("\n");
}

void print_binary(unsigned int num)
{
unsigned int mask=0x80;   //mask = [0000 0000 1000 0000]

while(mask > 0)
   {
   if((num & mask) == 0 )
         printf("0");
   else
         printf("1");
  mask = mask >> 1 ;  // Right Shift
   }
}
