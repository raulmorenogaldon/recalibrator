
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


char CIGAR1[15] = "45M3D1I20X";
char CIGAR1_l = 10;

char CIGAR1_elems[3] = {45, 3, 1};
char CIGAR1_type[3] = {'M', 'D', 'I'};
char CIGAR1_elem_l = 3;

char CIGAR2[10] = "25Z3D1D0M";
char CIGAR2_l = 9;

char CIGAR2_elems[2] = {3, 1};
char CIGAR2_type[2] = {'D', 'D'};
char CIGAR2_elem_l = 2;

#endif /* CHECK_AUX_LIBRARY_H_ */
