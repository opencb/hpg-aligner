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

#include "aux_vector.h"

/**
 * Initializes vector with initial values.
 */
ERROR_CODE
initialize_vector(uint32_t *vector, const size_t size, const uint32_t value)
{
	int i;

	if(size == 0)
		return INVALID_INPUT_SIZE_0;

	for(i = 0; i < size; i++)
	{
		vector[i] = value;
	}

	return NO_ERROR;
}

/**
 * Return vector of integers with initial values.
 */
ERROR_CODE
new_vector_uint32(const size_t size, const uint32_t value, uint32_t **out_vector)
{
	int i;
	uint32_t *vector;

	if(size == 0)
	{
		*out_vector = NULL;
		return INVALID_INPUT_SIZE_0;
	}

	vector = (uint32_t *)malloc(size * sizeof(uint32_t));

	for(i = 0; i < size; i++)
		vector[i] = value;

	//printf("Created new vector uint32 with %d positions, total size %lu bytes\n", (int)size, size * sizeof(uint32_t));

	*out_vector = vector;

	return NO_ERROR;
}

/**
 * Return vector of double with initial values.
 */
ERROR_CODE
new_vector_double(const size_t size, const double value, double **out_vector)
{
	int i;
	double *vector;

	if(size == 0)
	{
		*out_vector = NULL;
		return INVALID_INPUT_SIZE_0;
	}

	vector = (double *)malloc(size * sizeof(double));

	for(i = 0; i < size; i++)
		vector[i] = value;

	//printf("Created new vector double with %d positions, total size %lu bytes\n", (int)size, size * sizeof(double));

	*out_vector = vector;

	return NO_ERROR;
}

/**
 * Return vector of double with initial values.
 */
ERROR_CODE
max_value(double *vector, size_t size, double *max)
{
	double maximum = -DBL_MAX;
	double aux;
	int i;

	for(i = size-1; i >= 0; i--)
	{
		aux = vector[i];
		if(aux > maximum)
			maximum = aux;
	}

	*max = maximum;

	return NO_ERROR;
}

/**
 * Get component index of vector max value.
 */
ERROR_CODE
max_index(double *vector, size_t size, uint16_t *max_i)
{
	double max = -DBL_MAX;
	int index = -1;
	double aux;
	int i;

	for(i = size-1; i >= 0; i--)
	{
		aux = vector[i];
		if(aux > max)
		{
			max = aux;
			index = i;
		}
	}

	*max_i = index;

	return NO_ERROR;
}
