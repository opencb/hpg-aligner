#ifndef SW1_H
#define SW1_H

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <immintrin.h>
//#include "smith_waterman.h"
#include "aligners/sw/smith_waterman.h"
#include "aligners/sw/macros.h"
#include "aligners/sw/sse.h"
#include "aligners/sw/emboss.h"

#ifdef __AVX__
#define SIMD_DEPTH 8
#define SIMD_ALIGN 32
#else
#define SIMD_DEPTH 4
#define SIMD_ALIGN 16
#endif // __AVX__

//====================================================================================
// Smith-Waterman structures and functions (SIMD version)
//====================================================================================

/**
 * @brief Input structure for Smith-Waterman algorithm.
 *
 * Input structure based-on SIMD for the Smith-Waterman algorithm,
 * basically, it contains the pointers to the sequences to align,
 * the number of pointers depends on the SIMD 'depth', e.g.,
 * for SSE, SIMD depth is 4, for AVX, SIMD depth is 8.
 */
typedef struct sw_simd_input {
     unsigned int depth; /**< SIMD depth, e.g., 4 for SSE, 8 for AVX. */
     char** seq_p;       /**< Pointers to the target sequences. */
     char** ref_p;       /**< Pointers to the reference sequences. */
     unsigned int* seq_len_p; /**< Pointer to the target sequences lengths. */
     unsigned int* ref_len_p; /**< Pointer to the reference sequences lengths. */
} sw_simd_input_t;

//------------------------------------------------------------------------------------

/**
 * @brief Constructor for the @a sw_simd_input_t structure.
 * @param depth SIMD depth, e.g. 4 for SSE, 8 for AVX
 * @return Pointer to the new structure.
 *
 * @a sw_simd_input_t constructor that allocates memory for
 * the sequences and lengths pointers, the number of pointers
 * depends on the depth.
 */
sw_simd_input_t* sw_simd_input_new(unsigned int depth);

//------------------------------------------------------------------------------------

/**
 * @brief Destructor for the @a sw_simd_input_t structure.
 * @param input_p[out] pointer to the structure to free
 *
 * @a sw_simd_input_t destructor that frees the memory previously
 * allocated by the constructor @a sw_simd_input_new.
 */
void sw_simd_input_free(sw_simd_input_t* input_p);

//------------------------------------------------------------------------------------

/**
 * @brief Adds a new sequence pair (target/refernce) to be align.
 * @param seq_p pointer to the target sequence
 * @param seq_len length of the target sequence
 * @param ref_p pointer to the reference sequence
 * @param ref_len lenth of the reference sequence
 * @param index position inside the @a sw_simd_input_t structure, ranging from 0 to depth
 * @param input_p[out] pointer to the structure where the pair is inserted
 *
 * Adds the target/reference sequences pair to the @a input_p structure at the 
 * position indicated by @a index.
 */
void sw_simd_input_add(char* seq_p, unsigned int seq_len, char* ref_p, unsigned int ref_len, 
		       unsigned int index, sw_simd_input_t* input_p);

//------------------------------------------------------------------------------------

/**
 * @brief Displays a @a sw_simd_input_t structure.
 * @param depth position until displaying
 * @param input_p pointer to the structure storing the sequences
 *
 * Displays the sequences to align from position 0 to
 * position @a depth, in pairs: target/reference.
 */
void sw_simd_input_display(unsigned int depth, sw_simd_input_t* input_p);

//------------------------------------------------------------------------------------

/**
 * @brief Output structure for Smith-Waterman algorithm.
 *
 * Output structure based-on SIMD for the Smith-Waterman algorithm,
 * basically, it contains the pointers to the aligned sequences,
 * the number of pointers depends on the SIMD 'depth', e.g.,
 * for SSE, SIMD depth is 4, for AVX, SIMD depth is 8.
 */
typedef struct sw_simd_output {
  unsigned int depth;  /**< SIMD depth, e.g., 4 for SSE, 8 for AVX. */
  char** mapped_seq_p; /**< Pointers to the aligned target sequences. */
  char** mapped_ref_p; /**< Pointers to the aligned reference sequences. */
  unsigned int* mapped_len_p; /**< Pointers to the aligned target sequences lengths. */
  unsigned int* start_p;
  unsigned int* start_seq_p;
  float* score_p;      /**< Pointers to the resulting scores. */
  float* norm_score_p; /**< Pointers to the resulting normalized scores (0..1). */
} sw_simd_output_t;

//------------------------------------------------------------------------------------

/**
 * @brief Constructor for the @a sw_simd_output_t structure.
 * @param depth SIMD depth, e.g. 4 for SSE, 8 for AVX
 * @return Pointer to the new structure.
 *
 * @a sw_simd_output_t constructor that allocates memory for
 * the aligned sequences and lengths pointers, the number of pointers
 * depends on the depth.
 */
sw_simd_output_t* sw_simd_output_new(unsigned int depth);

//------------------------------------------------------------------------------------

/**
 * @brief Destructor for the @a sw_simd_output_t structure.
 * @param output_p[out] pointer to the structure to free
 *
 * @a sw_simd_output_t destructor that frees the memory previously
 * allocated by the constructor @a sw_simd_output_new.
 */
void sw_simd_output_free(sw_simd_output_t* output_p);

//------------------------------------------------------------------------------------

/**
 * @brief Displays a @a sw_simd_output_t structure.
 * @param depth position until displaying
 * @param output_p pointer to the structure storing the aligned sequences
 *
 * Displays the aligned sequences from position 0 to
 * position @a depth, in pairs: target/reference.
 */
void sw_simd_output_display(unsigned int depth, sw_simd_output_t* output_p);

//====================================================================================
// Smith-Waterman structures and functions (SIMD version)
//====================================================================================

typedef struct sw_simd_context {
     float gap_open;   /**< Penalty for gap openning. */
     float gap_extend; /**< Penalty for gap extending. */
#ifdef __AVX__
     __m256 zero_simd;       /**< Register set to zero (SIMD registers). */
     __m256 gap_open_simd;   /**< Penalty for gap openning (SIMD registers). */
     __m256 gap_extend_simd; /**< Penalty for gap extending (SIMD registers). */
#else
     __m128 zero_simd;       /**< Register set to zero (SIMD registers). */
     __m128 gap_open_simd;   /**< Penalty for gap openning (SIMD registers). */
     __m128 gap_extend_simd; /**< Penalty for gap extending (SIMD registers). */
#endif // __AVX__

     int x_size;   /**< x-length of the H score matrix. */
     int y_size;   /**< y-length of the H score matrix. */
     int max_size; /**< H score matrix size, i.e, x_size * y_size. */
          
     float substitution[2]; /**< Array for storing the match and mismatch penalties. */
     float matrix[128][128]; /**< Array for storing the match and mismatch penalties. */
     
     float *E; /**< E vector. */
     float *F; /**< F vector. */
     float *H; /**< H score matrix. */

     int *C; /**< Path direction pointer array, to traceback. */
     char *q_aux;
     char *r_aux;
     int aux_size;
     int H_size;
     int F_size; 

     int *compass_p; /**< Path direction pointer array, to traceback. */

     int seq_x_end[SIMD_DEPTH]; /**< x-positions of the maximum scores. */
     int ref_y_end[SIMD_DEPTH]; /**< y-positions of the maximum scores. */
     
     float *h_end[SIMD_DEPTH];  /**< pointers to the H cells with the maximum scores. */
     char *seq_end[SIMD_DEPTH]; /**< pointers to the chars of target sequences corresponding to the maximum scores. */
     char *ref_end[SIMD_DEPTH]; /**< pointers to the chars of reference sequences corresponding to the maximum scores. */
     int *compass_end[SIMD_DEPTH]; /**< pointers to the path direction pointer array corresponding to the maximum scores. */
     
     char *a_map; /**< temporary pointer to the aligned target sequence. */
     char *b_map; /**< temporary pointer to the aligned reference sequence. */    
} sw_simd_context_t;

//------------------------------------------------------------------------------------

/**
 * @brief Constructor for the @a sw_simd_context_t structure.
 * @param match penalty for match (e.g., 5.0)
 * @param mismatch penalty for mismatch (e.g., -4.0)
 * @param gap_open penalty for gap openning (e.g., 10.0)
 * @param gap_extend penalty for gap extending (e.g., 0.5)
 * @return Pointer to the new structure.
 *
 * @a sw_context_see_t constructor that allocates memory for
 * the Smith-Waterman context for the SIMD version. The SIMD context
 * consists of public parameters (match, mismatch, gap_open and gap_extend) and
 * private parameters (E, F, H,..).
 */
sw_simd_context_t* sw_simd_context_new(float match, float mismatch, float gap_open, float gap_extend);

//------------------------------------------------------------------------------------

/**
 * @brief Destructor for the @a sw_simd_context_t structure.
 * @param context_p[out] pointer to the structure to free
 *
 * @a sw_simd_context_t destructor that frees the memory previously
 * allocated by the constructor @a sw_context_see_t new.
 */
void sw_simd_context_free(sw_simd_context_t *context_p);

//------------------------------------------------------------------------------------

/**
 * @brief Updates the SIMD context structure.
 * @param x_size x-length of the score matrix
 * @param y_size y-length of the score matrix
 * @param context_p[out] pointer to the structure to update
 *
 * @a sw_simd_context_t destructor that frees the memory previously
 * allocated by the constructor @a sw_context_see_t new. This is
 * a private function, it should be called by the @a smithwaterman_see
 * function. Mainly, updates the E, F, H vectors sizes and allocate memory
 * for them.
 */
void sw_simd_context_update(int x_size, int y_size, sw_simd_context_t* context_p);

//------------------------------------------------------------------------------------

/**
 * @brief Performs the Smith-Waterman algorithm based-on SIMD instructions.
 * @param input_p pointer to the input sequences to align
 * @param[out] output_p pointer to the output aligned sequences
 * @param[out] context_p pointer to the SIMD context
 *
 * Based-on SIMD instrunctions, this function performs the Smith-Waterman algorithm,
 * it will perform 4 x Smith-Waterman in parallel using the XMM registers.
 */
void smith_waterman_simd(sw_simd_input_t* input_p, sw_simd_output_t* output_p, sw_simd_context_t* context_p);

//====================================================================================
// Smith-Waterman functions from EMBOSS package
//====================================================================================

#define LEFT 1
#define DOWN 2
#define U_FEPS 1.192e-6F         // 1.0F + E_FEPS != 1.0F
#define E_FPEQ(a,b,e) (((b - e) < a) && (a < (b + e)))

//-------------------------------------------------------------

float smith_waterman(char* seq_a, char* seq_b, float gapopen, float gapextend, char* m, char* n, int* start1, int* start2);

//------------------------------------------------------------------------------------
/*
float AlignPathCalcSW(const char *a, const char *b, int lena, int lenb,
                      float gapopen, float gapextend, float *path,
                      int *compass);

//------------------------------------------------------------------------------------

void AlignWalkSWMatrix(const float *path, const int *compass,
		       float gapopen, float gapextend,
		       const char*  a, const char* b,
		       char* m, char* n,
		       int lena, int lenb,
		       int *start1, int *start2);
//------------------------------------------------------------------------------------

void revstr(char* str);
*/
#endif // SW1_H

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
