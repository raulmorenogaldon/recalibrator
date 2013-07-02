/*
 * aux_bam.h
 *
 *  Created on: Jun 25, 2013
 *      Author: rmoreno
 */

#ifndef AUX_BAM_H_
#define AUX_BAM_H_

#include <bam_file.h>
#include <alignment.h>
#include <recal_common.h>

extern bam_header_t* create_empty_bam_header(unsigned int num_chroms);

extern void compare_bams_qual(const char* bamPath0, const char* bamPath1, int cycles);

#endif /* AUX_BAM_H_ */
