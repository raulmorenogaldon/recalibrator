#include <timestats.h>
#include <time.h>
#include <sys/timeb.h>
#include <stdint.h>

pthread_mutex_t time_mutex = PTHREAD_MUTEX_INITIALIZER;

p_timestats TIME_GLOBAL_STATS;

typedef struct time_slot {
	uint64_t sec;	//Internal use
	uint64_t nsec;
	uint64_t min_sec;	//Min time
	uint64_t min_nsec;
	uint64_t max_sec;	//Max time
	uint64_t max_nsec;
	uint64_t sum_sec;	//Total time
	uint64_t sum_nsec;
	uint64_t number;	//Number of times counted (for mean sum/number)
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
		
		slot->sec = 0;
		slot->nsec = 0;
		slot->min_sec = UINT64_MAX;
		slot->min_nsec = UINT64_MAX;
		slot->max_sec = 0;
		slot->max_nsec = 0;
		slot->sum_sec = 0;
		slot->sum_nsec = 0;
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

	if(!stats)
	{
		printf("Time - WARNING: Attempting to destroy NULL pointer time\n");
		return INVALID_INPUT_PARAMS_NULL;
	}

	time_stats_t *s = (time_stats_t *)(*stats);

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

	free(s->slots);
	free(s);
	*stats = NULL;

	return NO_ERROR;
}

/**
 * TIME OPERATIONS
 */
ERROR_CODE
time_init_slot(const unsigned int slot, p_timestats stats)
{
	time_stats_t *s = (time_stats_t *)stats;
	struct timespec ts;

	//Get time
	clock_gettime(CLOCK_REALTIME, &ts);

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
	s->slots[slot]->sec = (uint64_t)ts.tv_sec;
	s->slots[slot]->nsec = (uint64_t)ts.tv_nsec;
	pthread_mutex_unlock(&time_mutex);
}


ERROR_CODE
time_set_slot(const unsigned int slot, p_timestats stats)
{
	time_stats_t *s = (time_stats_t *)stats;
	struct timespec ts;
	uint64_t interval_sec, interval_nsec;
	
	clock_gettime(CLOCK_REALTIME, &ts);

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

	//Calc time intervals
	interval_sec = (uint64_t)ts.tv_sec - s->slots[slot]->sec;
	interval_nsec = (uint64_t)ts.tv_nsec - s->slots[slot]->nsec;

	pthread_mutex_lock(&time_mutex);

	if(s->slots[slot]->max_sec <= interval_sec && s->slots[slot]->max_nsec < interval_nsec)
	{
		s->slots[slot]->max_sec = interval_sec;
		s->slots[slot]->max_nsec = interval_nsec;
	}

	if(s->slots[slot]->min_sec >= interval_sec && s->slots[slot]->min_nsec > interval_nsec)
	{

		s->slots[slot]->min_sec = interval_sec;
		s->slots[slot]->min_nsec = interval_nsec;
	}
		
	s->slots[slot]->sum_sec += interval_sec;
	s->slots[slot]->sum_nsec += interval_nsec;
	s->slots[slot]->number++;
	
	pthread_mutex_unlock(&time_mutex);
}

ERROR_CODE
time_get_mean_slot(const unsigned int slot, const p_timestats stats, double *out_mean)
{
	time_stats_t *s = (time_stats_t *)stats;

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

	*out_mean = ( (double)s->slots[slot]->sum_sec + ((double) s->slots[slot]->sum_nsec / 1000000000.0) ) / (double)s->slots[slot]->number;

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

	*out_min = (double)s->slots[slot]->min_sec + ((double) s->slots[slot]->min_nsec / 1000000000.0);

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

	*out_max = (double)s->slots[slot]->max_sec + ((double) s->slots[slot]->max_nsec / 1000000000.0);

	return NO_ERROR;
}
