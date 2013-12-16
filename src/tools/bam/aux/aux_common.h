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

#ifndef AUX_COMMON_H_
#define AUX_COMMON_H_

#include <stdint.h>

#ifdef __MMX__
#include <mmintrin.h>
#endif

#ifdef __SSE__
#include <xmmintrin.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif

/**
 * NULL DEFINITION
 */

#ifndef NULL
	#ifdef __cplusplus
		#define NULL 0
	#else
		#define NULL (void *)0
	#endif
#endif

/**
 * EXTERN DEFINE
 */

#ifdef __cplusplus
	#define EXTERNC extern "C"
#else
	#define EXTERNC extern
#endif

/**
 * DATA LENGTHS
 */
typedef uint8_t U_QUALS;
typedef uint64_t U_BASES;
typedef uint32_t U_CYCLES;
typedef uint8_t U_DINUC;

/**
 * BOOLS
 */

#define BOOL unsigned char
#define TRUE 1
#define FALSE 0

/**
 * INLINING
 */

#ifndef __GNUC__
	#define __asm__ asm
	#define __inline__ inline
#endif
#define INLINE __inline__

/**
 * MEMORY ALIGNMENT
 */
#define MEM_ALIG_SIZE 16

/**
 * HOT AND COLD ATTRIBUTES
 */
#define __ATTR_HOT __attribute__((hot))

/**
 * ASSERTIONS
 */
#include <assert.h>
#define ASSERT(expr) assert(expr)

/**
 * ERROR CODES
 */

//#define ERROR_CODE unsigned char;
enum ERROR_C {
	NO_ERROR,

	INVALID_INPUT_PARAMS,
	INVALID_INPUT_PARAMS_NULL,
	INVALID_INPUT_PARAMS_0,
	INVALID_INPUT_PARAMS_NEGATIVE,
	INVALID_INPUT_SIZE_0,

	//Timestats specific
	INVALID_INPUT_SLOT,

	//Recalibrate specific
	INVALID_INPUT_QUAL,
	INVALID_INPUT_DINUC,
	ERROR_OBTAINING_SEQ,
	ERROR_OBTAINING_QUALS,
	INVALID_SEQ_LENGTH
};
typedef enum ERROR_C ERROR_CODE;

#endif /* AUX_COMMON_H_ */
