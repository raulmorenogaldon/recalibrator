
#ifndef CHECK_AUX_LIBRARY_H_
#define CHECK_AUX_LIBRARY_H_

#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include <time.h>
#include <recal_common.h>
#include <recal_config.h>
#include <aux_library.h>
#include <timestats.h>


//decompose_cigar
char CIGAR1[11] = "45M3D1I20X";
char CIGAR1_l = 10;

char CIGAR1_elems[4] = {45, 3, 1, 20};
char CIGAR1_type[4] = {'M', 'D', 'I', 'X'};
char CIGAR1_elem_l = 4;

char CIGAR2[10] = "25Z3D1D0M";
char CIGAR2_l = 9;

char CIGAR2_elems[3] = {25, 3, 1};
char CIGAR2_type[3] = {'Z', 'D', 'D'};
char CIGAR2_elem_l = 3;

//supress_indels
char seq[11] = "AGGTGCGCTT";
uint8_t seq_l = 10;
char seq_cigar[11] = "5=2I2=3D1M";
uint8_t seq_cigar_l = 10;
char seq_res[12] = "AGGTGCTXXXT";
uint8_t seq_res_l = 11;

#endif /* CHECK_AUX_LIBRARY_H_ */
