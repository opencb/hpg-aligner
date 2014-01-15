#ifndef AUX_MATH_H_
#define AUX_MATH_H_

#include "aux_library.h"


/***************************
 * MATH OPERATIONS
 **************************/

/**
 * Returns quality from probability.
 * \param p Error probability [0,1]
 */
EXTERNC double Qvalue(double p);

/**
 * Return probability from quality.
 * \param q Error quality. Default -> [P_SANGER_MIN,P_SANGER_MAX]
 */
EXTERNC double Pvalue(double q);

/**
 * Return Solexa quality from probability.
 * \param p Error probability [0,1]
 */
EXTERNC double Qsolexa(double p);

/**
 * Return probability from Solexa quality.
 * \param q Error quality [P_SOLEXA_MIN,P_SOLEXA_MAX]
 */
EXTERNC double Psolexa(double q);

/**
 * Return Sanger quality from probability.
 * \param p Error probability [0,1]
 */
EXTERNC double Qsanger(double p);

/**
 * Return probability from Sanger quality.
 * \param q Error quality [P_SANGER_MIN,P_SANGER_MAX]
 */
EXTERNC double Psanger(double q);

/**
 * Return evaluation of gaussian function.
 * g: a + b * exp(-pow(value-c,2)/(2*d*d))
 * \param value Value to evaluate.
 * \param a A constant
 * \param b B constant
 * \param c C constant
 * \param d D constant
 */
EXTERNC double gaussian_function(double value, double a, double b, double c, double d);

/**
 * Return approximation for log10(n!)
 * \param n Number to apply factorial and logarithm.
 */
EXTERNC double log10_gamma(uint64_t n);

#endif /* AUX_MATH_H_ */
