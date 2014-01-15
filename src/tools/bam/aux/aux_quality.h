
#ifndef AUX_QUALITY_H_
#define AUX_QUALITY_H_

#include "aux_library.h"

/***************************
 * QUALITY OPERATIONS
 **************************/

/**
 * Obtains estimated global quality from a vector of bases.
 * \param v_bases 	Vector containing number of bases for each quality.
 * 					Each component represents a count per each quality, beginning from start_quality value.
 * 					For i = 0 -> number of bases for quality start_quality.
 * 					For i = 1 -> number of bases for quality start_quality + 1 ...
 * \param count	Number of vector components.
 * \param start_quality Indicates Quality value for first component of v_bases.
 * \param estimated_Q Pointer to store result of estimated global quality.
 */
EXTERNC ERROR_CODE recal_get_estimated_Q(U_BASES *v_bases, size_t count, U_QUALS start_quality, double *estimated_Q);

/**
 * Obtains empirical quality from a count of bases and misses and a theorical quality.
 * \param miss Counts how many errors are in bases.
 * \param bases Counts how many bases.
 * \param initial_quality Theorical quality.
 * \param emp_qual Pointer to store result of empirical quality.
 */
EXTERNC ERROR_CODE recal_get_empirical_Q(double miss, U_BASES bases, double initial_quality, double *emp_qual);

/**
 * Obtains how much one empirical quality approximates one theorical quality, expressed by logarithm.
 * \param Qemp Empirical quality to be evaluated.
 * \param Qreported Theorical quality to compare.
 * \param log Pointer to store logarithm of comparation result.
 */
EXTERNC ERROR_CODE log10_Qemp_Reported(double Qemp, double Qreported, double *log);

/**
 * Obtains how much one empirical quality approximates one theorical quality expressed by a count of observations and errors.
 * Result is expressed by a logarithm.
 * \param Qemp Empirical quality to be evaluated.
 * \param obs Count of observations.
 * \param err Count of errors in this observations.
 * \param log Pointer to store logarithm of comparation result.
 */
EXTERNC ERROR_CODE log10_Qemp_likelihood(double Qemp, U_BASES obs, U_BASES err, double *log);

#endif /* AUX_QUALITY_H_ */
