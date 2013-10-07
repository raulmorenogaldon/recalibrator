#include <timestats.h>
#include <time.h>
#include <sys/timeb.h>
#include <stdint.h>
#include <float.h>

pthread_mutex_t time_mutex = PTHREAD_MUTEX_INITIALIZER;

p_timestats TIME_GLOBAL_STATS;

typedef struct time_slot {
	//Accumulators
	double dmin_sec;
	double dmax_sec;
	double dsum_sec;

	//Internal use
	uint64_t sec;
	uint64_t nsec;

	uint64_t number;	//Number of times counted (for mean sum/number)
} time_slot_t;

typedef struct time_stats {
	FILE *f_output;
	unsigned int num_slots;
	time_slot_t **slots;
} time_stats_t;

/**
 * PRIVATE FUNCTIONS
 */
ERROR_CODE time_aux_add_time(const unsigned int slot, p_timestats stats, const double time);

ERROR_CODE time_new_stats(const unsigned int num_slots, p_timestats *out_timestats)
{
	time_stats_t *stats;
	time_slot_t *slot;
	int i;
	
	stats = (time_stats_t *) malloc(sizeof(time_stats_t));
	
	stats->f_output = NULL;
	stats->num_slots = num_slots;
	stats->slots = (time_slot_t **)malloc(num_slots * sizeof(time_slot_t *));
	
	//Initialize slots
	for(i = 0; i < num_slots; i++)
	{
		slot = (time_slot_t *) malloc(sizeof(time_slot_t));
		
		slot->dmin_sec = DBL_MAX;
		slot->dmax_sec = 0.0;
		slot->dsum_sec = 0.0;
		slot->sec = 0;
		slot->nsec = 0;
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

	//If output file
	if(s->f_output)
		fclose(s->f_output);

	free(s);
	*stats = NULL;

	return NO_ERROR;
}

/**
 * Time statistics output to file
 */
ERROR_CODE
time_set_output_file(const char *name, p_timestats *stats)
{
	if(!stats)
	{
		printf("Time - WARNING: Attempting to access NULL pointer time\n");
		return INVALID_INPUT_PARAMS_NULL;
	}

	time_stats_t *s = (time_stats_t *)(*stats);

	if(!s)
	{
		printf("Time - WARNING: Attempting to access NULL pointer time\n");
		return INVALID_INPUT_PARAMS_NULL;
	}

	if(s->f_output)
	{
		fclose(s->f_output);
	}

	s->f_output = fopen(name, "w");

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

	pthread_mutex_lock(&time_mutex);

	if(slot >= s->num_slots)
	{
		pthread_mutex_unlock(&time_mutex);

		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return INVALID_INPUT_SLOT;
	}

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
	double interval;
	
	clock_gettime(CLOCK_REALTIME, &ts);

	if(!s)
	{
		printf("Time - WARNING: Attempting to set slot from NULL pointer time\n");
		return INVALID_INPUT_PARAMS_NULL;
	}

	pthread_mutex_lock(&time_mutex);
	if(slot >= s->num_slots)
	{
		pthread_mutex_unlock(&time_mutex);

		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return INVALID_INPUT_SLOT;
	}
	pthread_mutex_unlock(&time_mutex);

	//Calc time intervals
	interval_sec = (uint64_t)ts.tv_sec - s->slots[slot]->sec;
	if(ts.tv_nsec < s->slots[slot]->nsec)
	{
		interval_nsec = s->slots[slot]->nsec - (uint64_t)ts.tv_nsec;
		interval_nsec = 1000000000.0 - interval_nsec;
		interval_sec--;
	}
	else
	{
		interval_nsec = (uint64_t)ts.tv_nsec - s->slots[slot]->nsec;
	}
	interval = (double)interval_sec + ((double)interval_nsec / 1000000000.0);

	time_aux_add_time(slot, stats, interval);
}

ERROR_CODE
time_add_time_slot(const unsigned int slot, p_timestats stats, const double time)
{
	time_stats_t *s = (time_stats_t *)stats;

	if(!s)
	{
		printf("Time - WARNING: Attempting to add time to NULL pointer time\n");
		return INVALID_INPUT_PARAMS_NULL;
	}

	pthread_mutex_lock(&time_mutex);
	if(slot >= s->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return INVALID_INPUT_SLOT;
	}
	pthread_mutex_unlock(&time_mutex);

	if(time < 0)
	{
		printf("Time: Trying to add negative time = %lf\n", time);
		return INVALID_INPUT_PARAMS_NEGATIVE;
	}

	time_aux_add_time(slot, stats, time);
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

	pthread_mutex_lock(&time_mutex);
	if(slot >= s->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return INVALID_INPUT_SLOT;
	}
	pthread_mutex_unlock(&time_mutex);

	*out_mean = s->slots[slot]->dsum_sec / (double)s->slots[slot]->number;

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

	pthread_mutex_lock(&time_mutex);
	if(slot >= s->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return INVALID_INPUT_SLOT;
	}
	pthread_mutex_unlock(&time_mutex);

	*out_min = s->slots[slot]->dmin_sec;

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

	pthread_mutex_lock(&time_mutex);
	if(slot >= s->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return INVALID_INPUT_SLOT;
	}
	pthread_mutex_unlock(&time_mutex);

	*out_max = s->slots[slot]->dmax_sec;

	return NO_ERROR;
}

/**
 * PRIVATE FUNCTIONS
 */
ERROR_CODE
time_aux_add_time(const unsigned int slot, p_timestats stats, const double time)
{
	time_stats_t *s = (time_stats_t *)stats;

	pthread_mutex_lock(&time_mutex);

	if(s->slots[slot]->dmax_sec <= time)
	{
		s->slots[slot]->dmax_sec = time;
	}

	if(s->slots[slot]->dmin_sec >= time)
	{
		s->slots[slot]->dmin_sec = time;
	}

	s->slots[slot]->number++;
	s->slots[slot]->dsum_sec += time;

	//File output
	if(s->f_output)
	{
		//< SLOT TIME(s) >
		fprintf(s->f_output, "%d %lf\n", slot, time);
	}

	pthread_mutex_unlock(&time_mutex);
}
