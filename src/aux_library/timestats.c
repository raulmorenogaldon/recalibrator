#include <timestats.h>

pthread_mutex_t time_mutex = PTHREAD_MUTEX_INITIALIZER;

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

void *time_new_stats(int num_slots)
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

	return (void *)stats;
}

void time_init_slot(int slot, clock_t initial_time, void *stats)
{
	time_stats_t *s = (time_stats_t *)stats;

	if(s == NULL)
	{
		printf("Time - WARNING: Attempting to initialize slot from NULL pointer time\n");
		return;
	}
	if(slot > s->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return;
	}
	
	pthread_mutex_lock(&time_mutex);	
	s->slots[slot]->aux_time = initial_time;
	pthread_mutex_unlock(&time_mutex);
}


void time_set_slot(int slot, clock_t end_time, void *stats)
{
	double time;
	time_stats_t *s = (time_stats_t *)stats;
	
	if(s == NULL)
	{
		printf("Time - WARNING: Attempting to set slot from NULL pointer time\n");
		return;
	}
	if(slot > s->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return;
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

double time_get_mean_slot(int slot, void *stats)
{
	time_stats_t *s = (time_stats_t *)s;

	if(s == NULL)
	{
		printf("Time - WARNING: Attempting to get slot mean from NULL pointer time\n");
		return -1;
	}

	return (double) (s->slots[slot]->sum / s->slots[slot]->number);
}

double time_get_min_slot(int slot, void *stats)
{
	time_stats_t *s = (time_stats_t *)stats;

	if(s == NULL)
	{
		printf("Time - WARNING: Attempting to get slot min from NULL pointer time\n");
		return -1;
	}

	return s->slots[slot]->min;
}

double time_get_max_slot(int slot, void *stats)
{
	time_stats_t *s = (time_stats_t *)stats;

	if(s == NULL)
	{
		printf("Time - WARNING: Attempting to get slot max from NULL pointer time\n");
		return -1;
	}

	return s->slots[slot]->max;
}

void time_destroy_stats(void *stats)
{
	int i;
	time_stats_t *s = (time_stats_t *)stats;
	
	if(s == NULL)
	{
		printf("Time - WARNING: Attempting to destroy NULL pointer time\n");
		return -1;
	}

	//Destroy slots
	for(i = 0; i < s->num_slots; i++)
	{	
		free(s->slots[i]);
	}
	
	free(s);
}
