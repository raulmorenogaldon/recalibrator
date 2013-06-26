/*
 * aux_bam.h
 *
 *  Created on: Jun 25, 2013
 *      Author: rmoreno
 */

#ifndef AUX_BAM_H_
#define AUX_BAM_H_

#include "common.h"
#include "alignment.h"
#include "bam_file.h"
#include "bam.h"

void compare_bams_qual(const char* bamPath0, const char* bamPath1, int cycles);

bam_header_t* create_empty_bam_header(int num_chroms);

#endif /* AUX_BAM_H_ */
