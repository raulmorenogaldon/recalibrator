#ifndef POOL_H
#define POOL_H

#include <stdio.h>

#include "bioformats/bam-sam/alignment.h"
#include "bioformats/bam-sam/bam_file.h"


typedef struct bam_pool {
	unsigned int num_batchs;
	size_t batch_size;
	unsigned int pivot;
	unsigned int valid;
	void* pool;
}bam_pool_t;

bam_pool_t* pool_new(int num_batchs, size_t batch_size);
void pool_free(bam_pool_t *pool);


int pool_has_free_batchs(bam_pool_t *pool);
void pool_add_batch(bam_batch_t *batch, bam_pool_t *pool);
bam_batch_t *pool_get_batch(bam_pool_t *pool);
void pool_next_batch(bam_pool_t *pool);
#endif
