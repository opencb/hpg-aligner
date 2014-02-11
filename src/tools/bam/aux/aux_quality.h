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

#ifndef AUX_QUALITY_H_
#define AUX_QUALITY_H_

#include "aux_library.h"

/***************************
 * QUALITY OPERATIONS
 **************************/

/**
 * Obtains estimated global quality from a vector of bases.
 * \param[in] v_bases 	Vector containing number of bases for each quality.
 * 					Each component represents a count per each quality, beginning from start_quality value.
 * 					For i = 0 -> number of bases for quality start_quality.
 * 					For i = 1 -> number of bases for quality start_quality + 1 ...
 * \param[in] count	Number of vector components.
 * \param[in] start_quality Indicates Quality value for first component of v_bases.
 * \param[out] estimated_Q Pointer to store result of estimated global quality.
 */
EXTERNC ERROR_CODE recal_get_estimated_Q(U_BASES *v_bases, size_t count, U_QUALS start_quality, double *estimated_Q);

/**
 * Obtains empirical quality from a count of bases and misses and a theorical quality.
 * \param[in] miss Counts how many errors are in bases.
 * \param[in] bases Counts how many bases.
 * \param[in] initial_quality Theorical quality.
 * \param[out] emp_qual Pointer to store result of empirical quality.
 */
EXTERNC ERROR_CODE recal_get_empirical_Q(double miss, U_BASES bases, double initial_quality, double *emp_qual);

/**
 * Obtains how much one empirical quality approximates one theorical quality, expressed by logarithm.
 * \param[in] Qemp Empirical quality to be evaluated.
 * \param[in] Qreported Theorical quality to compare.
 * \param[out] log Pointer to store logarithm of comparation result.
 */
EXTERNC ERROR_CODE log10_Qemp_Reported(double Qemp, double Qreported, double *log);

/**
 * Obtains how much one empirical quality approximates one theorical quality expressed by a count of observations and errors.
 * Result is expressed by a logarithm.
 * \param[in] Qemp Empirical quality to be evaluated.
 * \param[in] obs Count of observations.
 * \param[in] err Count of errors in this observations.
 * \param[out] log Pointer to store logarithm of comparation result.
 */
EXTERNC ERROR_CODE log10_Qemp_likelihood(double Qemp, U_BASES obs, U_BASES err, double *log);

#endif /* AUX_QUALITY_H_ */
