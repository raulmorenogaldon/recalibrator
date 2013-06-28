#ifndef TIMESTATS_H
#define TIMESTATS_H

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

#include "common.h"

/**
 * Global time statistics for use
 */
static void *TIME_GLOBAL_STATS;

/**
 * Time statistics structure creation
 */
extern void *time_new_stats(int num_slots);

/**
 * Time statistics structure delete
 */
extern void time_destroy_stats(void *stats);

/**
 * TIME OPERATIONS
 */
extern void time_init_slot(int slot, clock_t initial_time, void *stats);
extern void time_set_slot(int slot, clock_t end_time, void *stats);
extern double time_get_mean_slot(int slot, void *stats);
extern double time_get_min_slot(int slot, void *stats);
extern double time_get_max_slot(int slot, void *stats);

#endif
