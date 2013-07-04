#ifndef TIMESTATS_H
#define TIMESTATS_H

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

/**
 * Time structure
 */
typedef void *p_timestats;

/**
 * Global time statistics for use
 */
extern p_timestats TIME_GLOBAL_STATS;

/**
 * Time statistics structure creation
 */
extern p_timestats time_new_stats(unsigned int num_slots);

/**
 * Time statistics structure delete
 */
extern void time_destroy_stats(p_timestats stats);

/**
 * TIME OPERATIONS
 */
extern void time_init_slot(unsigned int slot, clock_t initial_time, p_timestats stats);
extern void time_set_slot(unsigned int slot, clock_t end_time, p_timestats stats);
extern double time_get_mean_slot(unsigned int slot, p_timestats stats);
extern double time_get_min_slot(unsigned int slot, p_timestats stats);
extern double time_get_max_slot(unsigned int slot, p_timestats stats);

#endif
