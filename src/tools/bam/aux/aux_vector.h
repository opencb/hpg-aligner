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

#ifndef AUX_VECTOR_H_
#define AUX_VECTOR_H_

#include "aux_library.h"

/***************************
 * VECTOR OPERATIONS
 **************************/

/**
 * Initializes vector with initial values.
 */
EXTERNC ERROR_CODE initialize_vector(uint32_t *vector, const size_t size, const uint32_t value);

/**
 * Return vector of integers with initial values.
 */
EXTERNC ERROR_CODE new_vector_uint32(const size_t size, const uint32_t value, uint32_t **out_vector);

/**
 * Return vector of double with initial values.
 */
EXTERNC ERROR_CODE new_vector_double(const size_t size, const double value, double **out_vector);

/**
 * Get max value of vector.
 * \param[in] vector Vector to evaluate.
 * \param[in] size Size of vector.
 * \param[out] max Pointer to store maximum value.
 */
EXTERNC ERROR_CODE max_value(double *vector, size_t size, double *max);

/**
 * Get component index of vector max value.
 * \param[in] vector Vector to evaluate.
 * \param[in] size Size of vector.
 * \param[out] max_i Pointer to store maximum value component index.
 */
EXTERNC ERROR_CODE max_index(double *vector, size_t size, uint16_t *max_i);

#endif /* AUX_VECTOR_H_ */
