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

#ifndef AUX_NUCLEOTIDE_H_
#define AUX_NUCLEOTIDE_H_

#include "aux_library.h"
#include "x86intrin.h"

/***************************
 * NUCLEOTIDE OPERATIONS
 **************************/

/**
 * Compare two sequences and obtain missmatches. 
 * \param[in] ref_seq First sequence.
 * \param[in] bam_seq Second sequence. 
 * \param[in] bam_seq_l Length of input sequences. 
 * \param[out] comp_seq Output sequence to store individual missmatch value. 0 means equal, 1 means differ.
 * \param[out] miss_count Number of missmatches in comparation. Optional, set to NULL in case.
 */
static inline ERROR_CODE nucleotide_compare(const char *ref_seq, const char *bam_seq, size_t bam_seq_l, char *comp_res, uint32_t *miss_count) __ATTR_HOT __ATTR_INLINE;

/**
 * Compare two sequences and obtain missmatches and missmatches qualities summatory.
 * \param[in] ref_seq First sequence.
 * \param[in] bam_seq Second sequence.
 * \param[in] bam_qual Sequence qualities.
 * \param[in] bam_seq_l Length of input sequences.
 * \param[out] comp_seq Output sequence to store individual missmatch value. 0 means equal, 1 means differ.
 * \param[out] out_miss_count Number of missmatches in comparation. Optional, set to NULL in case.
 * \param[out] out_sum_quals Summatory of missmatch qualities.
 */
static inline ERROR_CODE nucleotide_miss_qual_sum(const char *ref_seq, const char *bam_seq, const char *bam_qual, size_t bam_seq_l, char *comp_res, uint32_t *out_miss_count, uint32_t *out_sum_quals) __ATTR_HOT __ATTR_INLINE;



/**
 * INLINE DEFINITIONS
 */

/**
 * Compare two sequences and obtain missmatches.
 */
static inline ERROR_CODE
nucleotide_compare(const char *ref_seq, const char *bam_seq, size_t bam_seq_l, char *comp_res, uint32_t *miss_count)
{
	int i;
	uint32_t count;
	uint32_t misses;

	assert(ref_seq);
	assert(bam_seq);
	assert(bam_seq_l > 0);
	assert(comp_res);

	//SSE
#ifdef __SSE2__	 //SSE2 block
	__m128i v_ref, v_seq, v_count;

	//Comparation
	__m128i v_comp;

	//Constants
	const __m128i v_zero = _mm_set1_epi8(0);
	const __m128i v_one = _mm_set1_epi8(1);
	const __m128i v_one16 = _mm_set1_epi16(1);

	char *ref_aux;
	char *bam_aux;
	char *comp_aux;
	size_t aux_l;

	//Allocate buffers
	aux_l = bam_seq_l + (16 - (bam_seq_l % 16));
	ref_aux = (char *) _mm_malloc(aux_l * sizeof(char), MEM_ALIG_SSE_SIZE);
	bam_aux = (char *) _mm_malloc(aux_l * sizeof(char), MEM_ALIG_SSE_SIZE);
	comp_aux = (char *) _mm_malloc(aux_l * sizeof(char), MEM_ALIG_SSE_SIZE);

	//Init buffers
	memcpy(ref_aux, ref_seq, bam_seq_l);
	memcpy(bam_aux, bam_seq, bam_seq_l);
	memset(ref_aux + bam_seq_l, 0, aux_l - bam_seq_l);	//Set last elements to match
	memset(bam_aux + bam_seq_l, 0, aux_l - bam_seq_l);

	//Init count
	v_count = _mm_set1_epi16(0);

	//Iterates nucleotides in this read
	misses = 0;
	for(i = 0; i < bam_seq_l; i += 16)
	{
		//Pack sequences
		v_ref = _mm_load_si128((__m128i const *)(ref_aux + i));
		v_seq = _mm_load_si128((__m128i const *)(bam_aux + i));

		//Compare sequences
		v_comp = _mm_cmpeq_epi8(v_ref, v_seq);	//0xFF Equals, 0x00 Diff

		//Change result
		v_comp = _mm_add_epi8(v_comp, v_one);	//Now equals = 0 and diff = 1

		//Sum number of differences
		//Accumulate 8 bit values to two 64 bit partial sum
		v_count = _mm_add_epi64(v_count, _mm_sad_epu8(v_comp, v_zero));

		//Store comparation values
		_mm_store_si128((__m128i*)(comp_aux + i), v_comp);
	}

	//Horizontal add of two 64 bit partial sums
	misses = _mm_cvtsi128_si32(v_count) + _mm_cvtsi128_si32(_mm_srli_si128(v_count, 8));

	//Set result
	memcpy(comp_res, comp_aux, bam_seq_l);

	//Set miss count
	if(miss_count)
	{
		*miss_count = misses;
	}

	//Free
	_mm_free(ref_aux);
	_mm_free(bam_aux);
	_mm_free(comp_aux);

#else //Sequential Block

	//Iterates nucleotides in this read
	count = 0;
	for(i = 0; i < bam_seq_l; i++)
	{
		//0 Equals, 1 Diff
		if(ref_seq[i] != bam_seq[i])
		{
			comp_res[i] = 1;
			count++;
		}
		else
		{
			comp_res[i] = 0;
		}

	}

	//Count misses
	if(miss_count)
	{
		*miss_count = count;
	}

#endif	//End SSE if

	return NO_ERROR;
}

/**
 * Compare two sequences and obtain missmatches and missmatches qualities summatory. 
 */
static inline ERROR_CODE
nucleotide_miss_qual_sum(const char *ref_seq, const char *bam_seq, const char *bam_qual, size_t bam_seq_l, char *comp_res, uint32_t *out_miss_count, uint32_t *out_sum_quals)
{
	int i, z;
	uint32_t misses;
	uint32_t sum;

	assert(ref_seq);
	assert(bam_seq);
	assert(bam_qual);
	assert(bam_seq_l > 0);
	assert(comp_res);

	//SSE
#ifdef __SSE2__	 //SSE2 block
	//Constants
	const __m128i v_zero = _mm_set1_epi8(0);
	const __m128i v_one = _mm_set1_epi8(1);
	const __m128i v_one16 = _mm_set1_epi16(1);

	//Nucleotides
	__m128i v_ref, v_seq, v_qual;
	__m128i vl_qual, vh_qual;

	//Comparations
	__m128i v_comp, vl_comp, vh_comp;

	//Sums
	__m128i v_sum, v_count;

	char *ref_aux;
	char *bam_aux;
	char *qual_aux;
	char *comp_aux;
	size_t aux_l;

	//Allocate buffers
	aux_l = bam_seq_l + (16 - (bam_seq_l % 16));
	ref_aux = (char *) _mm_malloc(aux_l * sizeof(char), MEM_ALIG_SSE_SIZE);
	bam_aux = (char *) _mm_malloc(aux_l * sizeof(char), MEM_ALIG_SSE_SIZE);
	qual_aux = (char *) _mm_malloc(aux_l * sizeof(char), MEM_ALIG_SSE_SIZE);
	comp_aux = (char *) _mm_malloc(aux_l * sizeof(char), MEM_ALIG_SSE_SIZE);

	//Init buffers
	memcpy(ref_aux, ref_seq, bam_seq_l);
	memcpy(bam_aux, bam_seq, bam_seq_l);
	memcpy(qual_aux, bam_qual, bam_seq_l);
	memset(ref_aux + bam_seq_l, 0, aux_l - bam_seq_l);	//Set last elements to match
	memset(bam_aux + bam_seq_l, 0, aux_l - bam_seq_l);
	memset(qual_aux + bam_seq_l, 0, aux_l - bam_seq_l);

	//Init count
	v_sum = _mm_set1_epi32(0);
	v_count = _mm_set1_epi32(0);

	//Iterates nucleotides in this read
	for(i = 0; i < bam_seq_l; i += 16)
	{
		//Pack sequences
		v_ref = _mm_load_si128((__m128i const *)(ref_aux + i));
		v_seq = _mm_load_si128((__m128i const *)(bam_aux + i));
		v_qual = _mm_load_si128((__m128i const *)(qual_aux + i));

		//Compare sequences
		v_comp = _mm_cmpeq_epi8(v_ref, v_seq);	//0xFF Equals, 0x00 Diff

		//Transform result
		v_comp = _mm_add_epi8(v_comp, v_one);	//Now equals = 0 and diff = 1

		//Sum number of differences
		v_count = _mm_add_epi64(v_count, _mm_sad_epu8(v_comp, v_zero));

		//Sum different qualities
		{
			//Convert diffs to 16 bits integer
			vl_comp = _mm_unpacklo_epi8(v_comp, v_zero);
			vh_comp = _mm_unpackhi_epi8(v_comp, v_zero);

			//Convert qualities to 16 bits integer
			vl_qual = _mm_unpacklo_epi8(v_qual, v_zero);
			vh_qual = _mm_unpackhi_epi8(v_qual, v_zero);

			//Accumulate 16 bit values to 32 bit partial sum (only different)
			v_sum = _mm_add_epi32(v_sum, _mm_madd_epi16(vl_qual, vl_comp));
			v_sum = _mm_add_epi32(v_sum, _mm_madd_epi16(vh_qual, vh_comp));
		}

		//Store comparation values
		_mm_store_si128((__m128i*)(comp_aux + i), v_comp);
	}

	//Horizontal add of four 32 bit partial sums
#ifdef __SSSE3__
	v_sum = _mm_hadd_epi32(_mm_hadd_epi32(v_sum, v_zero), v_zero);
#else
	v_sum = _mm_add_epi32(v_sum, _mm_srli_si128(v_sum, 8));
	v_sum = _mm_add_epi32(v_sum, _mm_srli_si128(v_sum, 4));
#endif
	sum = _mm_cvtsi128_si32(v_sum);

	v_count = _mm_add_epi32(v_count, _mm_srli_si128(v_count, 8));
	misses = _mm_cvtsi128_si32(v_count);

	//Set result
	memcpy(comp_res, comp_aux, bam_seq_l);

	//Free
	_mm_free(ref_aux);
	_mm_free(bam_aux);
	_mm_free(qual_aux);
	_mm_free(comp_aux);

#else //Sequential Block

	//Iterates nucleotides in this read
	misses = 0;
	sum = 0;
	for(i = 0; i < bam_seq_l; i++)
	{
		//0 Equals, 1 Diff
		if(ref_seq[i] != bam_seq[i])
		{
			comp_res[i] = 1;
			sum += bam_qual[i];
			misses++;
		}
		else
		{
			comp_res[i] = 0;
		}

	}

#endif	//End SSE if

	//Set miss count
	if(out_miss_count)
	{
		*out_miss_count = misses;
	}

	//Sum diff qualities
	if(out_sum_quals)
	{
		//Set output
		*out_sum_quals = sum;
	}

	return NO_ERROR;
}



#endif /* AUX_NUCLEOTIDE_H_ */
