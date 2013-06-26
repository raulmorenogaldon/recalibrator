#ifndef TIMESTATS_H
#define TIMESTATS_H

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

#include "common.h"

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


time_stats_t *time_global_stats;

time_stats_t *time_new_stats(int num_slots);

void time_init_slot(int slot, clock_t initial_time, time_stats_t *stats);
void time_set_slot(int slot, clock_t end_time, time_stats_t *stats);
double time_get_mean_slot(int slot, time_stats_t *stats);
double time_get_min_slot(int slot, time_stats_t *stats);
double time_get_max_slot(int slot, time_stats_t *stats);

void time_destroy_stats(time_stats_t *stats);

#endif
