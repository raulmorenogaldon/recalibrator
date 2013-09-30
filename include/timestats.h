#ifndef TIMESTATS_H
#define TIMESTATS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include <recal_common.h>


/**
 * Time structure
 */
typedef void * p_timestats;

/**
 * Global time statistics for use
 */
EXTERNC p_timestats TIME_GLOBAL_STATS;

/**
 * Time statistics structure creation
 */
EXTERNC ERROR_CODE time_new_stats(const unsigned int num_slots, p_timestats *out_timestats);

/**
 * Time statistics structure delete
 */
EXTERNC ERROR_CODE time_destroy_stats(p_timestats *stats);

/**
 * TIME OPERATIONS
 */
EXTERNC ERROR_CODE time_init_slot(const unsigned int slot, p_timestats stats);
EXTERNC ERROR_CODE time_set_slot(const unsigned int slot, p_timestats stats);
EXTERNC ERROR_CODE time_get_mean_slot(const unsigned int slot, const p_timestats stats, double *out_mean);
EXTERNC ERROR_CODE time_get_min_slot(const unsigned int slot, const p_timestats stats, double *out_min);
EXTERNC ERROR_CODE time_get_max_slot(const unsigned int slot, const p_timestats stats, double *out_max);

#endif
