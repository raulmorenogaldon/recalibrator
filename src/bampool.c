#include "bampool.h"


bam_pool_t* pool_new(int num_batchs, size_t batch_size)
{
	bam_pool_t *pool;
	
	printf("Creating pool with %lu bytes\n",batch_size*num_batchs);
	printf("Size of pool: %d\n", num_batchs);
	printf("Size of batch in pool: %lu bytes\n", batch_size);
	
	pool = malloc(sizeof(bam_pool_t));
	
	pool->num_batchs = num_batchs;
	pool->batch_size = batch_size;
	pool->pivot = 0;
	pool->valid = 0;
	
	pool->pool = (void *)malloc(num_batchs * batch_size);
	
	return pool;
}


void pool_free(bam_pool_t *pool)
{
	printf("Deleting pool\n");
	free(pool->pool);
	free(pool);
}



int pool_has_free_batchs(bam_pool_t *pool)
{
	return pool->num_batchs - pool->valid;
}

void pool_add_batch(bam_batch_t *batch, bam_pool_t *pool)
{
	int pos;
	void *pointer;
	int size = 0;
	int sizeaux = 0;
	int i;
	bam_batch_t *newbatch;
	
	//Pool has free batchs?
	if(!pool_has_free_batchs(pool))
		return;
		
	printf("\nAdding batch to pool...\n");
		
	//Get free batch
	pos = (pool->pivot + pool->valid) % pool->num_batchs;
	
	//Copy batch struct to pool
	pointer = (void *)(pool->pool + pos * pool->batch_size);
	newbatch = (bam_batch_t *)pointer;
	memcpy(pointer, batch, sizeof(bam_batch_t));
	printf("Size of batch struct: %lu\n", sizeof(bam_batch_t));
	size += sizeof(bam_batch_t);
	
	//Copy bam1_t list to pool
	pointer = (void *)(pointer + sizeof(bam_batch_t));
	memcpy(pointer, batch->alignments_p, sizeof(bam1_t*) * batch->allocated_alignments);
	newbatch->alignments_p = (bam1_t**)pointer;
	printf("Size of bam1_t list: %lu\n", sizeof(bam1_t*) * batch->allocated_alignments);
	size += sizeof(bam1_t*) * batch->allocated_alignments;
	
	//Copy bam1_t structs to pool
	pointer = (void *)(pointer + sizeof(bam1_t*) * batch->allocated_alignments);
	for(i = 0; i < batch->num_alignments; i++)
	{
		memcpy(pointer, batch->alignments_p[i], sizeof(bam1_t));
		newbatch->alignments_p[i] = (bam1_t*)pointer;
		pointer = (void *)(pointer + sizeof(bam1_t));
		size += sizeof(bam1_t);
	}
	printf("Size of bam1_t struct: %lu\n", sizeof(bam1_t));
	printf("Size of all bam1_t structs: %lu\n", sizeof(bam1_t) * batch->num_alignments);
	
	//Copy data to pool
	for(i = 0; i < batch->num_alignments; i++)
	{
		memcpy(pointer, batch->alignments_p[i]->data, sizeof(uint8_t) * batch->alignments_p[i]->m_data);
		newbatch->alignments_p[i]->data = (uint8_t*)pointer;
		pointer = (void *)(pointer + sizeof(uint8_t) * batch->alignments_p[i]->data_len);
		sizeaux += sizeof(uint8_t) * batch->alignments_p[i]->data_len;
	}
	printf("Size of bam1_t data: %d\n", sizeaux);
	printf("Total size: %d\n\n", size + sizeaux);

	pool->valid++;
}

bam_batch_t *pool_get_batch(bam_pool_t *pool)
{
	return (bam_batch_t *)(pool->pool + pool->pivot*pool->batch_size);
}

void pool_next_batch(bam_pool_t *pool)
{
	if(pool->valid > 1)
	{
		pool->valid--;
		pool->pivot = (pool->pivot + 1) % pool->num_batchs;
	}
}
