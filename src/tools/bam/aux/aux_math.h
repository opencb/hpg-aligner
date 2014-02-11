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

#ifndef AUX_MATH_H_
#define AUX_MATH_H_

#include "aux_library.h"


/***************************
 * MATH OPERATIONS
 **************************/

/**
 * Returns quality from probability.
 * \param[in] p Error probability [0,1]
 */
EXTERNC double Qvalue(double p);

/**
 * Return probability from quality.
 * \param[in] q Error quality. Default -> [P_SANGER_MIN,P_SANGER_MAX]
 */
EXTERNC double Pvalue(double q);

/**
 * Return Solexa quality from probability.
 * \param[in] p Error probability [0,1]
 */
EXTERNC double Qsolexa(double p);

/**
 * Return probability from Solexa quality.
 * \param[in] q Error quality [P_SOLEXA_MIN,P_SOLEXA_MAX]
 */
EXTERNC double Psolexa(double q);

/**
 * Return Sanger quality from probability.
 * \param[in] p Error probability [0,1]
 */
EXTERNC double Qsanger(double p);

/**
 * Return probability from Sanger quality.
 * \param[in] q Error quality [P_SANGER_MIN,P_SANGER_MAX]
 */
EXTERNC double Psanger(double q);

/**
 * Return evaluation of gaussian function.
 * g: a + b * exp(-pow(value-c,2)/(2*d*d))
 * \param[in] value Value to evaluate.
 * \param[in] a A constant
 * \param[in] b B constant
 * \param[in] c C constant
 * \param[in] d D constant
 */
EXTERNC double gaussian_function(double value, double a, double b, double c, double d);

/**
 * Return approximation for log10(n!)
 * \param[in] n Number to apply factorial and logarithm.
 */
EXTERNC double log10_gamma(uint64_t n);

#endif /* AUX_MATH_H_ */
