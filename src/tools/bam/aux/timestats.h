/**
* Copyright (C) 2013 Raúl Moreno Galdón
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TIMESTATS_H
#define TIMESTATS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>

#ifndef EXTERNC
	#ifdef __cplusplus
		#define EXTERNC extern "C"
	#else
		#define EXTERNC extern
	#endif
#endif


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
EXTERNC int time_new_stats(const unsigned int num_slots, p_timestats *out_timestats);

/**
 * Time statistics structure delete
 */
EXTERNC int time_destroy_stats(p_timestats *stats);

/**
 * Time statistics output to file
 */
EXTERNC int time_set_output_file(const char *name, p_timestats stats);

/**
 * TIME OPERATIONS
 */
EXTERNC int time_init_slot(const unsigned int slot, p_timestats stats);
EXTERNC int time_set_slot(const unsigned int slot, p_timestats stats);
EXTERNC int time_add_time_slot(const unsigned int slot, p_timestats stats, const double time);
EXTERNC int time_get_mean_slot(const unsigned int slot, const p_timestats stats, double *out_mean);
EXTERNC int time_get_min_slot(const unsigned int slot, const p_timestats stats, double *out_min);
EXTERNC int time_get_max_slot(const unsigned int slot, const p_timestats stats, double *out_max);

#endif
