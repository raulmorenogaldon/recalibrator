/*
 * prueba.c
 * 
 */


#include <stdio.h>
#include <stdlib.h>

#include "bioformats/bam-sam/alignment.h"
#include "bioformats/bam-sam/bam_file.h"
#include "common-libs/containers/list.h"
#include "bamaux.h"
#include "recal.h"

void print_bam(char* bamPath)
{
	bam_file_t *bamF;
	bam1_t* bamAlig;
	bam_batch_t* bamBatch;
	alignment_t* aligAlig;
	
	printf("Abriendo BAM en \"%s\" ...\n", bamPath);
	bamF = bam_fopen(bamPath);
	printf("BAM abierto!...\n");
	
	printf("Nombre: \t%s\n", bamF->filename);
	printf("N alig: \t%d\n", bamF->num_alignments);
	printf("Modo:   \t%s\n", bamF->mode);
	
	bamBatch= bam_batch_new(1, SINGLE_CHROM_BATCH);
	bam_fread_max_size(bamBatch, 1, 1, bamF);
	
	//Obtengo el primer alignment del bam
	bamAlig = bamBatch->alignments_p[0];
	aligAlig = alignment_new_by_bam(bamAlig, 1);
	alignment_print(aligAlig);
	
	printf("Cerrando BAM...\n");
	bam_batch_free(bamBatch, 1);
	bam_fclose(bamF);
	
	printf("BAM cerrado.\n");
}

void crear_sub_bam(char* bamPath, int num_chrom, int num_alig, char* newBamPath)
{
	bam_header_t* newBamFileHeader;
	bam_file_t *bamF, *newBamFile;
	bam_batch_t* bamBatch;
	alignment_t* aligAlig;
	int i;
	
	printf("Abriendo BAM en \"%s\" ...\n", bamPath);
	bamF = bam_fopen(bamPath);
	printf("BAM abierto!...\n");
	
	printf("Nombre: \t%s\n", bamF->filename);
	printf("N alig: \t%d\n", bamF->num_alignments);
	printf("Modo:   \t%s\n", bamF->mode);
	
	bamBatch = bam_batch_new(10000000, SINGLE_CHROM_BATCH);
	bam_fread_max_size(bamBatch, 10000000, 1, bamF);
	
	printf("Creando nuevo BAM con %d cromosomas y %d alignments\n", num_chrom, num_alig);
	
	//newBamFileHeader = bam_header_init();
	newBamFileHeader = create_empty_bam_header(num_chrom);
	
	newBamFile = bam_fopen_mode(newBamPath, newBamFileHeader, "w");
	
	bam_fwrite_header(newBamFileHeader, newBamFile);
	
	for(i = 0; i < num_alig; i++)
	{
		bam_fwrite(bamBatch->alignments_p[i], newBamFile);
		printf("i = %d\n",i);
		//aligAlig = alignment_new_by_bam(bamBatch->alignments_p[i], 1);
		//alignment_print(aligAlig);
		//getchar();	
	}
	
	
	printf("BAM creado!\n");	
	
	//bam_batch_free(bamBatch, 1);
	bam_fclose(bamF);
	bam_fclose(newBamFile);	
}

void comparar_fichero(char* bamPath0, char* bamPath1, char* bamPath2)
{
	bam_file_t* bamFile0;
	bam_file_t* bamFile1;
	bam_file_t* bamFile2;
	bam_batch_t* bamBatch0;
	bam_batch_t* bamBatch1;
	bam_batch_t* bamBatch2;
	bam1_t* bamAlig;
	alignment_t* aligAlig0;
	alignment_t* aligAlig1;
	alignment_t* aligAlig2;
	int diff1, diff2, diff3, i;
	
	printf("Abriendo BAM original en \"%s\" ...\n", bamPath0);
	printf("Abriendo BAM recalibrado en \"%s\" ...\n", bamPath1);
	printf("Abriendo BAM recalibrado en \"%s\" ...\n", bamPath2);
	bamFile0 = bam_fopen(bamPath0);
	bamFile1 = bam_fopen(bamPath1);
	bamFile2 = bam_fopen(bamPath2);
	printf("BAM abiertos!...\n");
	
	//Establece el campo num_alignments
	//bam_fcount(bamFile1);
	//bam_fcount(bamFile2);
	
	printf("Nombres: \t%s <===> %s\n", bamFile1->filename, bamFile2->filename);
	printf("N aligs: \t%d <===> %d\n", bamFile1->num_alignments, bamFile2->num_alignments);
	printf("Modos:   \t%s <===> %s\n", bamFile1->mode, bamFile2->mode);
	
	
	/*bamBatch0 = bam_batch_new(9000, MULTIPLE_CHROM_BATCH);
	bamBatch1 = bam_batch_new(9000, MULTIPLE_CHROM_BATCH);
	bamBatch2 = bam_batch_new(9000, MULTIPLE_CHROM_BATCH);
	bam_fread_max_size(bamBatch0, 9000, 1, bamFile0);
	bam_fread_max_size(bamBatch1, 9000, 1, bamFile1);
	bam_fread_max_size(bamBatch2, 9000, 1, bamFile2);
	bam_batch_free(bamBatch0, 1);	
	bam_batch_free(bamBatch1, 1);
	bam_batch_free(bamBatch2, 1);*/
	bamBatch0 = bam_batch_new(1, MULTIPLE_CHROM_BATCH);
	bamBatch1 = bam_batch_new(1, MULTIPLE_CHROM_BATCH);
	bamBatch2 = bam_batch_new(1, MULTIPLE_CHROM_BATCH);
	bam_fread_max_size(bamBatch0, 1, 1, bamFile0);
	bam_fread_max_size(bamBatch1, 1, 1, bamFile1);
	bam_fread_max_size(bamBatch2, 1, 1, bamFile2);	
	
	//Obtengo el primer alignment del bam original
	bamAlig = bamBatch0->alignments_p[0];
	aligAlig0 = alignment_new_by_bam(bamAlig, 1);
	alignment_print(aligAlig0);
	
	//Obtengo el primer alignment del primer bam
	bamAlig = bamBatch1->alignments_p[0];
	aligAlig1 = alignment_new_by_bam(bamAlig, 1);
	alignment_print(aligAlig1);
	
	//Obtengo el primer alignment del segundo bam
	bamAlig = bamBatch2->alignments_p[0];
	aligAlig2 = alignment_new_by_bam(bamAlig, 1);
	alignment_print(aligAlig2);
	
	printf("\n\n---------------------------------------------------------\n");
	
	//Obtengo diferencia entre calidades
	printf("Diferencias: \nNuc\tQrecal1\tQref\tQrecal2\n");
	diff1=0;
	diff2=0;
	diff3=0;
	for(i=0; i < 100; i++)
	{
		printf("%c \t%d ", aligAlig1->sequence[i], aligAlig1->quality[i]);
		if(aligAlig1->quality[i] == aligAlig0->quality[i])
		{
			printf("====\t%d", aligAlig0->quality[i]);	
		}
		else
		{
			printf("\t%d", aligAlig0->quality[i]);
		}
		
		if(aligAlig0->quality[i] == aligAlig2->quality[i])
		{
			printf(" ====\t%d\n", aligAlig2->quality[i]);	
		}
		else
		{
			printf("\t%d\n", aligAlig2->quality[i]);
		}
		
		diff1 += abs(aligAlig1->quality[i] - aligAlig0->quality[i]);
		diff2 += abs(aligAlig2->quality[i] - aligAlig0->quality[i]);
		diff3 += abs(aligAlig2->quality[i] - aligAlig1->quality[i]);
	}
	printf("Diferencia entre calidades del BAM-1 total: %d\n", diff1);
	printf("Diferencia entre calidades del BAM-2 total: %d\n", diff2);
	printf("Diferencia entre calidades de BAM recalibrados: %d\n", diff3);
	
	printf("\n---------------------------------------------------------\n");
	printf("Cerrando BAM...\n");
	bam_fclose(bamFile0);
	bam_fclose(bamFile1);
	bam_fclose(bamFile2);
	bam_batch_free(bamBatch0, 1);
	bam_batch_free(bamBatch1, 1);
	bam_batch_free(bamBatch2, 1);
	
	printf("BAM cerrados.\n");
}

/*void probar_recalibrator_structs()
{
	recal_data_t* data;
	unsigned int qual, cycle;
	int miss;
	int i;
	
	//Creo un nuevo data
	data = recal_new_data();
	
	//Relleno con datos aleatorios
	for(i = 0; i < 10000; i++)
	{
		qual = rand() % 5 + 28; //28-32
		cycle = rand() % 5 + 1; //1-5
		miss = rand() % 2; //0-1
		add_cycle_base(data, qual, cycle, miss);
	}
	
	//Printeo el data
	fprint_data(data, "data_struct.txt");
	
	//Libero memoria
	destroy_data(data);
}*/

/*void probar_recalibrator(char* bam_path)
{
	recal_data_t* data;
	bam_file_t *bam_f;
	bam_file_t *bamF;
	
	//Creo un nuevo data
	data = recal_new_data();
	
	//Abro el bam
	printf("Abriendo BAM en \"%s\" ...\n", bam_path);
	bam_f = bam_fopen(bam_path);
	printf("BAM abierto!...\n");
	
	//Muestro datos
	//bam_fcount(bam_f);
	//printf("Nombre: \t%s\n", bam_f->filename);
	//printf("N alig: \t%d\n", bam_f->num_alignments);
	//printf("Modo:   \t%s\n", bam_f->mode);
	
	//Relleno el data
	get_data_from_bam(bam_f, data);
	
	//Printeo el data
	fprint_data(data, "recal_data.txt");
	
	//Libero memoria	
	printf("Cerrando BAM...\n");
	bam_fclose(bam_f);
	destroy_data(data);
	printf("BAM cerrado.\n");
}*/

void probar_bam_info(char* bam_path)
{
	bam_info_t* data;
	bam_file_t *bam_f;
	
	//Creo un nuevo data
	data = bam_new_info();
	
	//Abro el bam
	printf("Abriendo BAM en \"%s\" ...\n", bam_path);
	bam_f = bam_fopen(bam_path);
	printf("BAM abierto!...\n");
	
	//Relleno el data
	bam_get_info(bam_f, data);
	
	//Printeo el data
	printf("INFO: \n QUALS: min = %d, max = %d\n",data->min_qual, data->max_qual);
	printf("CYCLES: %d\n", data->num_cycles);
	
	//Libero memoria	
	printf("Cerrando BAM...\n");
	bam_fclose(bam_f);
	bam_destroy_info(data);
	printf("BAM cerrado.\n");
}

void probar_recal_struct()
{
	recal_info_t* data;
	int qual, cycle, dinuc;
	int miss;
	int i;
	
	//Creo un nuevo info
	data = recal_new_info();
	
	//Relleno con datos aleatorios
	for(i = 0; i < 1000000; i++)
	{
		qual = rand() % 5 + 28; //28-32
		cycle = rand() % 100; //0-99
		dinuc = rand() % 12; //0-11
		miss = rand() % 2; //0-1
		recal_add_base(data, qual, cycle, dinuc, miss);
	}
	
	//Printeo estadisticas
	printf("Tamaño contador: %lu (Bytes)\n", sizeof(unsigned int));
	printf("Tamaño vector calidad: %lu (Bytes)\n", (MAX_QUALITY - MIN_QUALITY) * sizeof(unsigned int));
	printf("Tamaño matriz calidad-ciclo: %lu (Bytes)\n", (MAX_QUALITY - MIN_QUALITY) * NUM_CYCLES * sizeof(unsigned int));
	printf("Tamaño matriz calidad-dinuc: %lu (Bytes)\n", (MAX_QUALITY - MIN_QUALITY) * NUM_DINUC * sizeof(unsigned int));
	
	//Printeo el data
	recal_fprint_info(data, "info_struct.txt");
	
	//Libero memoria
	recal_destroy_info(data);
}

void probar_recal(char* bam_path, char* ref_name, char* ref_path, char* info_path)
{
	recal_info_t* data;
	genome_t* ref;
	bam_file_t *bam_f;
	clock_t start, end;
	double cpu_time_used;
	  
    start = clock();
     
	//Creo un nuevo data
	data = recal_new_info();
	
	//Abro el bam
	printf("Abriendo BAM en \"%s\" ...\n", bam_path);
	bam_f = bam_fopen(bam_path);
	printf("BAM abierto!...\n");
	
	//Relleno la info del genoma referencia
	printf("Abriendo genoma referencia en \"%s%s\" ...\n", ref_path, ref_name);	
	ref = genome_new(ref_name, ref_path);
	printf("Genoma referencia abierto!...\n");
	
	//Relleno el data
	recal_get_data_from_bam(bam_f, ref, data);
	
	//Tiempo
	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	
	//Printeo el data
	recal_fprint_info(data, info_path);
	
	//Libero memoria	
	printf("Cerrando BAM...\n");
	bam_fclose(bam_f);
	recal_destroy_info(data);
	genome_free(ref);
	printf("BAM cerrado.\n");
	printf("Time elapsed: %.2f sec\n",cpu_time_used);
}

void probar_recal2(char* orig_bam_path, char* recal_bam_path, char* ref_name, char* ref_path, char* info_path)
{
	recal_info_t* data;
	
	//Time measures
	#ifdef D_TIME_DEBUG
		//Initialize stats
		time_new_stats(10);
	#endif
	
	//Creo un nuevo data
	data = recal_new_info();
	
	//Relleno el data
	recal_get_data_from_file(orig_bam_path, ref_name, ref_path, data);
	
	//Proceso deltas
	recal_calc_deltas(data);
	
	//Printeo el data
	recal_fprint_info(data, info_path);
	
	//Creo el bam
	recal_recalibrate_bam_file(orig_bam_path, data, recal_bam_path);
	
	#ifdef D_TIME_DEBUG
		//Print time stats
		printf("----------------------------\nTIME STATS: \n");
		
		printf("Time used to process BAM -> %.2f s\n", time_global_stats->slots[D_SLOT_GET_DATA_BAM]->min);
		printf("Time used to process one BATCH(mean) -> %.2f ms - min/max = %.2f/%.2f\n", 
			time_get_mean_slot(D_SLOT_GET_DATA_BATCH, time_global_stats)*1000.0, 
			time_get_min_slot(D_SLOT_GET_DATA_BATCH, time_global_stats)*1000.0,
			time_get_max_slot(D_SLOT_GET_DATA_BATCH, time_global_stats)*1000.0);
		printf("Time used to process one ALIGNMENT(mean) -> %.2f micros - min/max = %.2f/%.2f\n", 
			time_get_mean_slot(D_SLOT_GET_DATA_ALIG, time_global_stats)*1000000.0, 
			time_get_min_slot(D_SLOT_GET_DATA_ALIG, time_global_stats)*1000000.0,
			time_get_max_slot(D_SLOT_GET_DATA_ALIG, time_global_stats)*1000000.0);
		printf("Time used to read one BATCH(mean) -> %.2f ms - min/max = %.2f/%.2f\n", 
			time_get_mean_slot(D_SLOT_READ_BATCH, time_global_stats)*1000.0, 
			time_get_min_slot(D_SLOT_READ_BATCH, time_global_stats)*1000.0,
			time_get_max_slot(D_SLOT_READ_BATCH, time_global_stats)*1000.0);
		printf("Time used to process deltas -> %.2f ns\n", time_global_stats->slots[D_SLOT_CALC_DELTAS]->min*1000000000.0);
		printf("Time used to recalibrate -> %.2f s\n", time_global_stats->slots[D_SLOT_RECALIBRATE]->min);
			
		//Free memory from stats
		time_destroy_stats(time_global_stats);
	#endif
	
	//Libero memoria	
	recal_destroy_info(data);
}

void comprobar_memoria(char* bamPath)
{
	bam_file_t *bamF;
	bam1_t* bamAlig;
	bam_batch_t* bamBatch;
	alignment_t* aligAlig;
	alignment_t* aux_alig;
	int i;
	
	printf("Abriendo BAM en \"%s\" ...\n", bamPath);
	bamF = bam_fopen(bamPath);
	printf("BAM abierto!...\n");
	
	printf("Nombre: \t%s\n", bamF->filename);
	printf("N alig: \t%d\n", bamF->num_alignments);
	printf("Modo:   \t%s\n", bamF->mode);
	
	bamBatch= bam_batch_new(1, MULTIPLE_CHROM_BATCH);
	bam_fread_max_size(bamBatch, 1, 1, bamF);
	
	//Obtengo el primer alignment del bam
	bamAlig = bamBatch->alignments_p[0];
	
	for(i = 0; i < 100000; i++)
	{
		aux_alig = alignment_new_by_bam(bamAlig, 1);
		alignment_free(aux_alig);
	}
	
	printf("Cerrando BAM...\n");
	bam_batch_free(bamBatch, 1);
	bam_fclose(bamF);
	
	printf("BAM cerrado.\n");
}

int comprobar_data(recal_info_t* data1, recal_info_t* data2)
{
	int i,j;
	
	if(	data1->min_qual != data2->min_qual ||
		data1->num_quals != data2->num_quals ||
		data1->num_cycles != data2->num_cycles ||
		data1->num_dinuc != data2->num_dinuc)
	{
		return 1;
	}
	
	if(	data1->total_miss != data2->total_miss ||
		data1->total_bases != data2->total_bases ||
		data1->total_delta != data2->total_delta )
	{
		return 2;
	}
	
	for(i = 0; i < data1->num_quals; i++)
	{
		if(	data1->qual_miss[i] != data2->qual_miss[i] ||
			data1->qual_bases[i] != data2->qual_bases[i] ||
			data1->qual_delta[i] != data2->qual_delta[i] )
		{
			return 3;
		}
	}
	
	for(i = 0; i < data1->num_quals; i++)
	{
		for(j = 0; j < data1->num_cycles; j++)
		{
			if(	data1->qual_cycle_miss[i*data1->num_cycles + j] != data2->qual_cycle_miss[i*data1->num_cycles + j] ||
				data1->qual_cycle_bases[i*data1->num_cycles + j] != data2->qual_cycle_bases[i*data1->num_cycles + j] ||
				data1->qual_cycle_delta[i*data1->num_cycles + j] != data2->qual_cycle_delta[i*data1->num_cycles + j] )
			{
				return 4;
			}
		}
	}
	
	for(i = 0; i < data1->num_quals; i++)
	{
		for(j = 0; j < data1->num_dinuc; j++)
		{
			if(	data1->qual_dinuc_miss[i*data1->num_dinuc + j] != data2->qual_dinuc_miss[i*data1->num_dinuc + j] ||
				data1->qual_dinuc_bases[i*data1->num_dinuc + j] != data2->qual_dinuc_bases[i*data1->num_dinuc + j] ||
				data1->qual_dinuc_delta[i*data1->num_dinuc + j] != data2->qual_dinuc_delta[i*data1->num_dinuc + j] )
			{
				return 5;
			}
		}
	}
	
	return 0;
}

void probar_recal_paso1(char* orig_bam_path, char* ref_name, char* ref_path, char* info_path, char* info_bin_path)
{
	recal_info_t *data, *data2;
	int err;
	
	//Time measures
	#ifdef D_TIME_DEBUG
		//Initialize stats
		time_new_stats(10);
	#endif
	
	//Creo un nuevo data
	data = recal_new_info();
	data2 = recal_new_info();
	
	//Relleno el data
	recal_get_data_from_file(orig_bam_path, ref_name, ref_path, data);
	
	//Proceso deltas
	recal_calc_deltas(data);
	
	//Printeo el data
	recal_fprint_info(data, info_path);
	
	//Guardo los contadores
	recal_save_recal_info(data, info_bin_path);
	
	//Compruebo que el data del fichero es correcto
	recal_load_recal_info(info_bin_path, data2);
	
	err = comprobar_data(data, data2);
	if(err)
	{
		printf("ERROR: No coinciden los data = %d\n", err);
		recal_fprint_info(data2, "dataerr.txt");
	}
	else
	{
		printf("Guardado correcto\n");
	}
	
	#ifdef D_TIME_DEBUG
		//Print time stats
		printf("----------------------------\nTIME STATS: \n");
		
		printf("Time used to process BAM -> %.2f s\n", time_global_stats->slots[D_SLOT_GET_DATA_BAM]->min);
		printf("Time used to process one BATCH(mean) -> %.2f ms - min/max = %.2f/%.2f\n", 
			time_get_mean_slot(D_SLOT_GET_DATA_BATCH, time_global_stats)*1000.0, 
			time_get_min_slot(D_SLOT_GET_DATA_BATCH, time_global_stats)*1000.0,
			time_get_max_slot(D_SLOT_GET_DATA_BATCH, time_global_stats)*1000.0);
		printf("Time used to process one ALIGNMENT(mean) -> %.2f micros - min/max = %.2f/%.2f\n", 
			time_get_mean_slot(D_SLOT_GET_DATA_ALIG, time_global_stats)*1000000.0, 
			time_get_min_slot(D_SLOT_GET_DATA_ALIG, time_global_stats)*1000000.0,
			time_get_max_slot(D_SLOT_GET_DATA_ALIG, time_global_stats)*1000000.0);
		printf("Time used to read one BATCH(mean) -> %.2f ms - min/max = %.2f/%.2f\n", 
			time_get_mean_slot(D_SLOT_READ_BATCH, time_global_stats)*1000.0, 
			time_get_min_slot(D_SLOT_READ_BATCH, time_global_stats)*1000.0,
			time_get_max_slot(D_SLOT_READ_BATCH, time_global_stats)*1000.0);
		printf("Time used to process deltas -> %.2f ns\n", time_global_stats->slots[D_SLOT_CALC_DELTAS]->min*1000000000.0);
			
		//Free memory from stats
		time_destroy_stats(time_global_stats);
	#endif
	
	//Libero memoria	
	recal_destroy_info(data);
}

void probar_recal_paso2(char* orig_bam_path, char* recal_bam_path, char* ref_name, char* ref_path, char* info_bin_path)
{
	recal_info_t* data;
	
	//Time measures
	#ifdef D_TIME_DEBUG
		//Initialize stats
		time_new_stats(10);
	#endif
	
	//Creo un nuevo data
	data = recal_new_info();
	
	//Cargo el data
	recal_load_recal_info(info_bin_path, data);
	
	//Creo el bam
	recal_recalibrate_bam_file(orig_bam_path, data, recal_bam_path);
	
	#ifdef D_TIME_DEBUG
		//Print time stats
		printf("----------------------------\nTIME STATS: \n");
		printf("Time used to recalibrate -> %.2f s\n", time_global_stats->slots[D_SLOT_RECALIBRATE]->min);
			
		//Free memory from stats
		time_destroy_stats(time_global_stats);
	#endif
	
	//Libero memoria	
	recal_destroy_info(data);
}

//#define C_BAM_ERR08
//#define C_BAM_TEST
#define C_BAM_EXAMPLE


int main(int argc, char **argv)
{
	
	#ifdef C_BAM_ERR08
	char* path_bams = "../bams/recal_bams";
	char* bamPath = "../bams/recal_bams/ERR089819.bam";
	char* bamPathRecal = "../bams/recal_output/outputrecal.bam";
	char* bamPathGATK = "../bams/recal_output/outputgatk.bam";
	char* info_bin_name = "c_bam_err08.data";
	#endif
	
	#ifdef C_BAM_TEST
	char* path_bams = "../bams/recal_bams";
	char* bamPath = "../bams/recal_bams/testbam.bam";
	char* bamPathRecal = "../bams/recal_output/testoutputrecal.bam";
	char* bamPathGATK = "../bams/recal_output/testoutputgatk.bam";
	char* info_bin_name = "c_bam_test.data";
	#endif
	
	#ifdef C_BAM_EXAMPLE
	char* path_bams = "../bams/recal_bams";
	char* bamPath = "../bams/example_bams/exampleBAM.bam"; 
	char* bamPathRecal = "../bams/recal_output/exampleoutrecal.bam"; 
	char* bamPathGATK = "../bams/recal_output/exampleoutgatk.bam"; 
	char* info_bin_name = "c_bam_example.data";
	#endif
	
	char* ref_name = "/refpreprocesada.fa";
	char* info_name = "info_struct2.txt";
	
	init_log();
	#ifdef DEBUG	
	LOG_LEVEL(LOG_DEBUG_LEVEL);
	#endif
	
	//crear_sub_bam(bamPath, 1, 10, bamPathPrueba);
	
	//print_bam(bamPathPrueba);
	//print_bam(bamPath);
	
	//probar_recalibrator_structs();
	//probar_recalibrator(bamPath);
	//probar_bam_info(bamPath);
	//probar_recal_struct();
	//probar_recal(bamPath, ref_name, path_bams, info_name);
	
	
	//probar_recal2(bamPath, bamPathRecal, ref_name, path_bams, info_name);
	//comparar_fichero(bamPath, bamPathGATK, bamPathRecal);
	
	//probar_recal_paso1(bamPath, ref_name, path_bams, info_name, info_bin_name);
	//probar_recal_paso2(bamPath, bamPathRecal, ref_name, path_bams, info_bin_name);
	comparar_fichero(bamPath, bamPathGATK, bamPathRecal);
	
	//comprobar_memoria(bamPath);
	
	//recal_preprocess_fasta("../bams/example_bams/exampleFASTA.fasta", "../bams/example_bams/refexpreprocesada.fa");

	return 0;
}

