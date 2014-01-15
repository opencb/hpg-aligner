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
 * \param vector Vector to evaluate.
 * \param size Size of vector.
 * \param max Pointer to store maximum value.
 */
EXTERNC ERROR_CODE max_value(double *vector, size_t size, double *max);

/**
 * Get component index of vector max value.
 * \param vector Vector to evaluate.
 * \param size Size of vector.
 * \param max_i Pointer to store maximum value component index.
 */
EXTERNC ERROR_CODE max_index(double *vector, size_t size, uint16_t *max_i);

#endif /* AUX_VECTOR_H_ */
