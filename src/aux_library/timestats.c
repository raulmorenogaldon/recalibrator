#include "timestats.h"

pthread_mutex_t time_mutex = PTHREAD_MUTEX_INITIALIZER;

time_stats_t *time_new_stats(int num_slots)
{
	time_stats_t *stats;
	time_slot_t *slot;
	int i;
	
	stats = (time_stats_t *) malloc(sizeof(time_stats_t));
	
	stats->num_slots = num_slots;
	stats->slots = (time_slot_t **)malloc(num_slots * sizeof(time_slot_t *));
	
	//Initialize slots
	for(i = 0; i < num_slots; i++)
	{
		slot = (time_slot_t *) malloc(sizeof(time_slot_t));
		
		slot->min = 4000000000;
		slot->max = 0;
		slot->sum = 0;
		slot->number = 0;
		
		stats->slots[i] = slot;
	}
	
	time_global_stats = stats;
	
	return stats;
}

void time_init_slot(int slot, clock_t initial_time, time_stats_t *stats)
{
	if(stats == NULL)
	{
		return;
	}
	if(slot > stats->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", stats->num_slots);
		return;
	}
	
	pthread_mutex_lock(&time_mutex);	
	stats->slots[slot]->aux_time = initial_time;
	pthread_mutex_unlock(&time_mutex);
}


void time_set_slot(int slot, clock_t end_time, time_stats_t *stats)
{
	double time;
	
	if(stats == NULL)
	{
		return;
	}
	if(slot > stats->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", stats->num_slots);
		return;
	}
		
	pthread_mutex_lock(&time_mutex);
		
	time = ((double) (end_time - stats->slots[slot]->aux_time)) / CLOCKS_PER_SEC;
	
	if(stats->slots[slot]->max < time)
		stats->slots[slot]->max = time;
		
	if(stats->slots[slot]->min > time)
		stats->slots[slot]->min = time;
		
	stats->slots[slot]->sum += time;
	stats->slots[slot]->number++;
	
	pthread_mutex_unlock(&time_mutex);
}

double time_get_mean_slot(int slot, time_stats_t *stats)
{
	return (double) (stats->slots[slot]->sum / stats->slots[slot]->number);
}

double time_get_min_slot(int slot, time_stats_t *stats)
{
	return stats->slots[slot]->min;
}

double time_get_max_slot(int slot, time_stats_t *stats)
{
	return stats->slots[slot]->max;
}

void time_destroy_stats(time_stats_t *stats)
{
	int i;
	
	//Destroy slots
	for(i = 0; i < stats->num_slots; i++)
	{	
		free(stats->slots[i]);
	}
	
	free(stats);
}
