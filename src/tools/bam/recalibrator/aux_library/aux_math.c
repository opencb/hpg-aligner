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

#include "aux_math.h"

/**
 * Returns quality from probability.
 */
inline double
Qvalue(double P)
{
	#ifdef P_SOLEXA
		return Qsolexa(P);
	#else
		return Qsanger(P);
	#endif
}

/**
 * Return probability from quality.
 */
inline double
Pvalue(double Q)
{
	#ifdef P_SOLEXA
		return Psolexa(Q);
	#else
		return Psanger(Q);
	#endif
}

/**
 * Return Solexa quality from probability.
 */
inline double
Qsolexa(double p)
{
	double ret = -10.0 * log10(p/(1-p));

	if(ret > P_SOLEXA_MAX)
		ret = P_SOLEXA_MAX;
	if(ret < P_SOLEXA_MIN)
		ret = P_SOLEXA_MIN;

	return ret;
}

/**
 * Return probability from Solexa quality.
 */
inline double
Psolexa(double Q)
{
	double ret = pow(10.0, -Q/10.0) / (pow(10.0, -Q/10.0) + 1);

	if(ret > 1.0)
		ret = 1.0;
	if(ret < 0.0)
		ret = 0.0;

	return ret;
}

/**
 * Return Sanger quality from probability.
 */
inline double
Qsanger(double p)
{
	double ret = -10.0 * log10(p);

	if(ret > P_SANGER_MAX)
		ret = P_SANGER_MAX;
	if(ret < P_SANGER_MIN)
		ret = P_SANGER_MIN;

	return ret;
}

/**
 * Return probability from Solexa quality.
 */
inline double
Psanger(double Q)
{
	double ret = pow(10.0, -Q/10.0);

	if(ret > 1.0)
		ret = 1.0;
	if(ret < 0.0)
		ret = 0.0;

	return ret;
}

/**
 * Return evaluation of gaussian function.
 */
double
gaussian_function(double value, double a, double b, double c, double d)
{
	return a + b * exp(-pow(value-c,2)/(2*d*d));
}

/**
 * Return approximation for log10(n!)
 */
double
log10_gamma(uint64_t n)
{
	double ln = lgamma(n+1);

	double l10 = ln * log10(M_E);

	return l10;
}
