#include <timestats.h>

pthread_mutex_t time_mutex = PTHREAD_MUTEX_INITIALIZER;

p_timestats TIME_GLOBAL_STATS;

typedef struct time_slot {
	clock_t aux_time;	//Internal use
	double min;	//Min time
	double max;	//Max time
	double sum;	//Total time
	unsigned int number;	//Number of times counted (for mean sum/number)
} time_slot_t;

typedef struct time_stats {
	unsigned int num_slots;
	time_slot_t **slots;
} time_stats_t;

ERROR_CODE time_new_stats(const unsigned int num_slots, p_timestats *out_timestats)
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

	TIME_GLOBAL_STATS = stats;
	*out_timestats = stats;

	return NO_ERROR;
}

/**
 * Time statistics structure delete
 */
ERROR_CODE
time_destroy_stats(p_timestats *stats)
{
	int i;
	time_stats_t *s = (time_stats_t *)stats;

	if(!s)
	{
		printf("Time - WARNING: Attempting to destroy NULL pointer time\n");
		return INVALID_INPUT_PARAMS_NULL;
	}

	//Destroy slots
	for(i = 0; i < s->num_slots; i++)
	{
		free(s->slots[i]);
	}

	free(s);
	*stats = NULL;
}

/**
 * TIME OPERATIONS
 */
ERROR_CODE
time_init_slot(const unsigned int slot, const clock_t initial_time, p_timestats stats)
{
	time_stats_t *s = (time_stats_t *)stats;

	if(!s)
	{
		printf("Time - WARNING: Attempting to initialize slot from NULL pointer time\n");
		return INVALID_INPUT_PARAMS_NULL;
	}
	if(slot >= s->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return INVALID_INPUT_SLOT;
	}
	
	pthread_mutex_lock(&time_mutex);	
	s->slots[slot]->aux_time = initial_time;
	pthread_mutex_unlock(&time_mutex);
}


ERROR_CODE
time_set_slot(const unsigned int slot, const clock_t end_time, p_timestats stats)
{
	double time;
	time_stats_t *s = (time_stats_t *)stats;
	
	if(!s)
	{
		printf("Time - WARNING: Attempting to set slot from NULL pointer time\n");
		return INVALID_INPUT_PARAMS_NULL;
	}
	if(slot >= s->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return INVALID_INPUT_SLOT;
	}
		
	pthread_mutex_lock(&time_mutex);
		
	time = ((double) (end_time - s->slots[slot]->aux_time)) / CLOCKS_PER_SEC;
	
	if(s->slots[slot]->max < time)
		s->slots[slot]->max = time;
		
	if(s->slots[slot]->min > time)
		s->slots[slot]->min = time;
		
	s->slots[slot]->sum += time;
	s->slots[slot]->number++;
	
	pthread_mutex_unlock(&time_mutex);
}

ERROR_CODE
time_get_mean_slot(const unsigned int slot, const p_timestats stats, double *out_mean)
{
	time_stats_t *s = (time_stats_t *)s;

	if(!s)
	{
		printf("Time - WARNING: Attempting to get slot mean from NULL pointer time\n");
		return INVALID_INPUT_PARAMS_NULL;
	}

	if(slot >= s->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return INVALID_INPUT_SLOT;
	}

	*out_mean = (double) (s->slots[slot]->sum / s->slots[slot]->number);

	return NO_ERROR;
}

ERROR_CODE
time_get_min_slot(const unsigned int slot, const p_timestats stats, double *out_min)
{
	time_stats_t *s = (time_stats_t *)stats;

	if(!s)
	{
		printf("Time - WARNING: Attempting to get slot min from NULL pointer time\n");
		return INVALID_INPUT_PARAMS_NULL;
	}

	if(slot >= s->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return INVALID_INPUT_SLOT;
	}

	*out_min = s->slots[slot]->min;

	return NO_ERROR;
}

ERROR_CODE
time_get_max_slot(const unsigned int slot, const p_timestats stats, double *out_max)
{
	time_stats_t *s = (time_stats_t *)stats;

	if(!s)
	{
		printf("Time - WARNING: Attempting to get slot max from NULL pointer time\n");
		return INVALID_INPUT_PARAMS_NULL;
	}

	if(slot >= s->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return INVALID_INPUT_SLOT;
	}

	*out_max = s->slots[slot]->max;

	return NO_ERROR;
}
