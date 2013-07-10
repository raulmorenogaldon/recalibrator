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

EXTERNC ERROR_CODE create_empty_bam_header(const unsigned int num_chroms, bam_header_t *out_header);

EXTERNC ERROR_CODE compare_bams_qual(const char* bamPath0, const char* bamPath1, const int cycles);

#endif /* AUX_BAM_H_ */
