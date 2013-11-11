
#ifndef CHECK_RECAL_H_
#define CHECK_RECAL_H_

#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include <time.h>
#include <recal_common.h>
#include <recal_config.h>
#include <aux_library.h>
#include <bam_recal_library.h>
#include <timestats.h>

#include "bam_recal/recal_structs.h"

#define D_SIZE 33

//TEST DATA
const size_t DATA_SIZE = D_SIZE;
const double ESTIMATED_Q = 25.7403;
const U_QUALS MINIMUM_QUALITY = 6;

//Global
const double GLOBAL_MISS = 7147.44;
const U_BASES GLOBAL_OBS = 700787;
const U_BASES GLOBAL_EMPIRICAL = 20;
const double GLOBAL_DELTA = -5.740300000000001;

const U_QUALS quality[D_SIZE] = {
		6, 		7, 		8, 		9, 		10,
		11, 	12, 	13, 	14, 	15, 	16, 	17, 	18, 	19, 	20,
		21, 	22, 	23, 	24, 	25, 	26, 	27, 	28, 	29, 	30,
		31, 	32, 	33, 	34, 	35, 	36, 	37, 	38
};

const double qual_errors[D_SIZE] = {
		43.90, 	71.77, 	70.23, 	120.54, 71.51,
		93.09, 	297.29, 114.91, 7.68, 	30.17, 	29.38, 	41.03, 	43.66, 	36.58, 	217.82,
		97.83,	52.60,	41.44,	54.12,	209.29,	140.19,	42.39, 	271.76, 134.53, 279.56,
		176.62, 390.35, 319.97, 309.83, 625.71,	717.60, 632.62, 1361.47
};

const U_BASES qual_obs[D_SIZE] = {
		360, 	569, 	743, 	1385, 	1049, 											//6
		1646, 	5471, 	1870, 	597, 	1165, 	1149, 	1588, 	1719, 	2438, 	11675,	//11
		4196, 	3429, 	2553, 	3259, 	18129, 	7152, 	3951, 	22848, 	11205, 	30077, 	//21
		13786, 	37182, 	36969, 	37029, 	60503, 	88754, 	84749, 	201592					//31
};

const double qual_empirical[D_SIZE] = {
		18.0, 	16.0, 	17.0, 	15.0, 	17.0,
		17.0, 	14.0, 	16.0, 	20.0, 	19.0, 	19.0, 	19.0, 	19.0, 	19.0, 	18.0,
		18.0, 	19.0, 	19.0,	19.0, 	20.0, 	18.0, 	20.0, 	19.0, 	19.0, 	20.0,
		19.0, 	20.0, 	21.0,	21.0, 	20.0,	21.0, 	21.0, 	22.0
};

const double qual_deltas[D_SIZE] = {
		-2.0, 	-4.0, 	-3.0, 	-5.0, 	-3.0,
		-3.0, 	-6.0, 	-4.0, 	0.0, 	-1.0, 	-1.0, 	-1.0, 	-1.0, 	-1.0, 	-2.0,
		-2.0, 	-1.0, 	-1.0,	-1.0, 	0.0, 	-2.0, 	0.0, 	-1.0, 	-1.0, 	0.0,
		-1.0, 	0.0, 	1.0,	1.0, 	0.0,	1.0, 	1.0, 	2.0
};

#endif /* CHECK_RECAL_H_ */
