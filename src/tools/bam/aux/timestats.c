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

#include "timestats.h"
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
int time_aux_add_time(const unsigned int slot, p_timestats stats, const double time);

int
time_new_stats(const unsigned int num_slots, p_timestats *out_timestats)
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

	return 0;
}

/**
 * Time statistics structure delete
 */
int
time_destroy_stats(p_timestats *stats)
{
	int i;

	if(!stats)
	{
		printf("Time - WARNING: Attempting to destroy NULL pointer time\n");
		return 1;
	}

	time_stats_t *s = (time_stats_t *)(*stats);

	if(!s)
	{
		printf("Time - WARNING: Attempting to destroy NULL pointer time\n");
		return 1;
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

	return 0;
}

/**
 * Time statistics output to file
 */
int
time_set_output_file(const char *name, p_timestats stats)
{
	if(!stats)
	{
		printf("Time - WARNING: Attempting to access NULL pointer time\n");
		return 1;
	}

	time_stats_t *s = (time_stats_t *)stats;

	if(!s)
	{
		printf("Time - WARNING: Attempting to access NULL pointer time\n");
		return 1;
	}

	if(s->f_output)
	{
		fclose(s->f_output);
	}

	s->f_output = fopen(name, "a");

	return 0;
}

/**
 * TIME OPERATIONS
 */
int
time_init_slot(const unsigned int slot, p_timestats stats)
{
	time_stats_t *s = (time_stats_t *)stats;
	struct timespec ts;

	//Get time
	clock_gettime(CLOCK_REALTIME, &ts);

	if(!s)
	{
		printf("Time - WARNING: Attempting to initialize slot from NULL pointer time\n");
		return 1;
	}

	pthread_mutex_lock(&time_mutex);

	if(slot >= s->num_slots)
	{
		pthread_mutex_unlock(&time_mutex);

		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return 1;
	}

	s->slots[slot]->sec = (uint64_t)ts.tv_sec;
	s->slots[slot]->nsec = (uint64_t)ts.tv_nsec;
	pthread_mutex_unlock(&time_mutex);

	return 0;
}


int
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
		return 1;
	}

	pthread_mutex_lock(&time_mutex);
	if(slot >= s->num_slots)
	{
		pthread_mutex_unlock(&time_mutex);

		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return 1;
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

	return time_aux_add_time(slot, stats, interval);
}

int
time_add_time_slot(const unsigned int slot, p_timestats stats, const double time)
{
	time_stats_t *s = (time_stats_t *)stats;

	if(!s)
	{
		printf("Time - WARNING: Attempting to add time to NULL pointer time\n");
		return 1;
	}

	pthread_mutex_lock(&time_mutex);
	if(slot >= s->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return 1;
	}
	pthread_mutex_unlock(&time_mutex);

	if(time < 0)
	{
		printf("Time: Trying to add negative time = %lf\n", time);
		return 1;
	}

	return time_aux_add_time(slot, stats, time);
}

int
time_get_mean_slot(const unsigned int slot, const p_timestats stats, double *out_mean)
{
	time_stats_t *s = (time_stats_t *)stats;

	if(!s)
	{
		printf("Time - WARNING: Attempting to get slot mean from NULL pointer time\n");
		return 1;
	}

	pthread_mutex_lock(&time_mutex);
	if(slot >= s->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return 1;
	}
	pthread_mutex_unlock(&time_mutex);

	*out_mean = s->slots[slot]->dsum_sec / (double)s->slots[slot]->number;

	return 0;
}

int
time_get_min_slot(const unsigned int slot, const p_timestats stats, double *out_min)
{
	time_stats_t *s = (time_stats_t *)stats;

	if(!s)
	{
		printf("Time - WARNING: Attempting to get slot min from NULL pointer time\n");
		return 1;
	}

	pthread_mutex_lock(&time_mutex);
	if(slot >= s->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return 1;
	}
	pthread_mutex_unlock(&time_mutex);

	*out_min = s->slots[slot]->dmin_sec;

	return 0;
}

int
time_get_max_slot(const unsigned int slot, const p_timestats stats, double *out_max)
{
	time_stats_t *s = (time_stats_t *)stats;

	if(!s)
	{
		printf("Time - WARNING: Attempting to get slot max from NULL pointer time\n");
		return 1;
	}

	pthread_mutex_lock(&time_mutex);
	if(slot >= s->num_slots)
	{
		printf("Time: illegal slot, maximum = %d\n", s->num_slots);
		return 1;
	}
	pthread_mutex_unlock(&time_mutex);

	*out_max = s->slots[slot]->dmax_sec;

	return 0;
}

/**
 * PRIVATE FUNCTIONS
 */
int
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
		fprintf(s->f_output, "%d %.9f\n", slot, time);
	}

	pthread_mutex_unlock(&time_mutex);

	return 0;
}
