/*
 * aux_simd.h
 *
 *  Created on: Feb 14, 2014
 *      Author: rmoreno
 */

#ifndef AUX_SIMD_H_
#define AUX_SIMD_H_

#include "aux_library.h"

/**
 * SSE2 FUNCTIONS
 */

/**
 * \brief Compare two vectors of 8 bits elements and returns a vector with differences.
 * \param[in] v_1 First vector to compare.
 * \param[in] v_2 Second vector to compare.
 * \param[out] v_out_comp Output vector with differences for every element comparation. 0 means equals and 1 means missmatch.
 * \return Number of found differences
 */
static inline int
sse2_epi8_compare_and_count_diffs(__m128i v_1, __m128i v_2, __m128i *v_out_comp)
{
	char aux[16];
	__m128i v_comp, v_count;

	//Compare sequences
	v_comp = _mm_cmpeq_epi8(v_1, v_2);	//0xFF Equals, 0x00 Diff

	//Change result
	v_comp = _mm_add_epi8(v_comp, _mm_set1_epi8(1));	//Now equals = 0 and diff = 1

	//Sum differences
	v_count = _mm_add_epi8(v_comp, _mm_slli_si128(v_comp, 8));
	v_count = _mm_add_epi8(v_count, _mm_slli_si128(v_count, 4));
	v_count = _mm_add_epi8(v_count, _mm_slli_si128(v_count, 2));
	v_count = _mm_add_epi8(v_count, _mm_slli_si128(v_count, 1));	//First element is differences count

	//Output
	*v_out_comp = v_comp;

	return (int) ((char *)&v_count)[15];
}


#endif /* AUX_SIMD_H_ */
