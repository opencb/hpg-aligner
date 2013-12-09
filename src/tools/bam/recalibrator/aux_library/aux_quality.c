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

#include "aux_quality.h"

/***************************
 * QUALITY OPERATIONS
 **************************/

/**
 * Obtains estimated global quality from a vector of bases.
 */
ERROR_CODE
recal_get_estimated_Q(U_BASES *v_bases, size_t count, U_QUALS start_quality, double *estimated_Q)
{
	double err0;
	double total_sum;
	U_BASES total_bases;
	double aux;
	U_QUALS quality;
	int i;

	total_sum = 0.0;
	total_bases = 0.0;

	//Obtains global base count and error count
	for(i = count-1; i >= 0; i--)
	{
		quality = (U_QUALS)i + start_quality;	/* Get quality for this component*/
		err0 = (double)v_bases[i] * pow(10.0, (-((double)quality)*0.1));	/* Calculate error based in count of bases and what quality have */
		//Increment global counters
		total_sum += err0;
		total_bases = total_bases + v_bases[i];
	}

	//Calculate estimated global probability
	aux = total_sum / (double)total_bases;

	//Obtain quality from probability
	aux = Qvalue(aux);
	*estimated_Q = aux;	/* Set estimated global quality in result pointer */

	return NO_ERROR;
}

/**
 * Obtains empirical quality from a count of bases and misses and a theorical quality.
 */
ERROR_CODE
recal_get_empirical_Q(double miss, U_BASES bases, double initial_quality, double *emp_qual)
{
	U_BASES mismatches;
	U_BASES observations;
	double log10[91];
	double norm[91];
	double sum_norm;
	uint16_t mle;
	double qempprior;
	double qemplike;
	double max;
	int i;

	mismatches = (U_BASES)(miss + 0.5) + SMOOTH_CONSTANT;
	observations = bases + (SMOOTH_CONSTANT*2);

	//Get logarithms
	for(i = 90; i >= 0; i--)
	{
		log10_Qemp_Reported((double)i, initial_quality, &qempprior);
		log10_Qemp_likelihood((double)i, observations, mismatches, &qemplike);
		log10[i] =  qempprior + qemplike;
	}

	//Normalize
	sum_norm = 0.0;
	max_value(log10, 91, &max);
	for(i = 90; i >= 0; i--)
	{
		norm[i] = pow(10, log10[i] - max);
		sum_norm += norm[i];
	}
	for(i = 90; i >= 0; i--)
	{
		norm[i] = norm[i] / sum_norm;
	}

	//Get maximum log index
	max_index(norm, 91, &mle);

	*emp_qual = (double)mle;

	return NO_ERROR;
}

/**
 * Obtains how much one empirical quality approximates one theorical quality, expressed by logarithm.
 */
ERROR_CODE
log10_Qemp_Reported(double Qemp, double Qreported, double *log)
{
	U_QUALS difference;
	double gaussian;
	double local_log;

	difference = abs((int)(Qemp - Qreported));

	if(difference > 40)
		difference = 40;

	gaussian = gaussian_function(difference, 0.0, 0.9, 0.0, 0.5);

	local_log = log10(gaussian);

	if(isinf(local_log))
		local_log = -DBL_MAX;

	*log = local_log;

	return NO_ERROR;
}

/**
 * Obtains how much one empirical quality approximates one theorical quality expressed by a count of observations and errors.
 * Result is expressed by a logarithm.
 */
ERROR_CODE
log10_Qemp_likelihood(double Qemp, U_BASES obs, U_BASES err, double *log)
{
	double log10p;
	double log10OneMinusP;

	if(obs == 0)
		return 0.0;

	log10p = -Qemp*0.1;


	log10OneMinusP = log10(1 - pow( 10, log10p) );

	if(isnan(log10OneMinusP) || isinf(log10OneMinusP))
		log10OneMinusP = -DBL_MAX;

	double factobs = log10_gamma(obs);
	double facterr = log10_gamma(err);
	double factdiff = log10_gamma(obs - err);

	*log = factobs - facterr - factdiff
			+ (log10p * err)
			+ (log10OneMinusP * (obs - err));

	return NO_ERROR;
}
