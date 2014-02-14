#include "aux_nucleotide.h"

#ifndef _REENTRANT
static char *ref_aux;
static char *bam_aux;
static char *comp_aux;
static size_t aux_l = 0;
#endif

/**
 * Compare two sequences and obtain missmatches.
 */
inline ERROR_CODE
nucleotide_compare(const char *ref_seq, const char *bam_seq, size_t bam_seq_l, char *comp_res, uint32_t *miss_count)
{
	int i;
	uint32_t count;
	uint32_t misses;

	assert(ref_seq);
	assert(bam_seq);
	assert(bam_seq_l > 0);
	assert(comp_res);
	//assert(miss_count);

	//SSE
#ifdef __SSE2__	 //SSE2 block
	__m128i v_ref, v_seq, v_comp, v_sum;

#ifdef _REENTRANT	//Reentrant block
	char *ref_aux;
	char *bam_aux;
	char *comp_aux;
	size_t aux_l;

	aux_l = bam_seq_l + (16 - (bam_seq_l % 16));
	ref_aux = (char *) _mm_malloc(aux_l * sizeof(char), MEM_ALIG_SSE_SIZE);
	bam_aux = (char *) _mm_malloc(aux_l * sizeof(char), MEM_ALIG_SSE_SIZE);
	comp_aux = (char *) _mm_malloc(aux_l * sizeof(char), MEM_ALIG_SSE_SIZE);

#else	//No reentrant block

	//Too small buffers?
	if(bam_seq_l > aux_l)
	{
		if(ref_aux)
			_mm_free(ref_aux);

		if(bam_aux)
			_mm_free(bam_aux);

		if(comp_aux)
			_mm_free(comp_aux);

		aux_l = bam_seq_l + (16 - (bam_seq_l % 16));
		ref_aux = (char *) _mm_malloc(aux_l * sizeof(char), MEM_ALIG_SSE_SIZE);
		bam_aux = (char *) _mm_malloc(aux_l * sizeof(char), MEM_ALIG_SSE_SIZE);
		comp_aux = (char *) _mm_malloc(aux_l * sizeof(char), MEM_ALIG_SSE_SIZE);
	}

#endif	//Reentrant end

	//Init buffers
	memcpy(ref_aux, ref_seq, bam_seq_l);
	memcpy(bam_aux, bam_seq, bam_seq_l);
	memset(ref_aux + bam_seq_l, 0, aux_l - bam_seq_l);	//Set last elements to match
	memset(bam_aux + bam_seq_l, 0, aux_l - bam_seq_l);

	//Iterates nucleotides in this read
	misses = 0;
	for(i = 0; i < bam_seq_l; i += 16)
	{
		//Prefetching
		__builtin_prefetch (ref_aux + i + 16, 0, 2);		//Read, high temporal locality
		__builtin_prefetch (bam_aux + i + 16, 0, 2);		//Read, high temporal locality
		__builtin_prefetch (comp_aux + i + 16, 1, 2);	//Write, high temporal locality

		//Pack sequences
		v_ref = _mm_load_si128(ref_aux + i);
		v_seq = _mm_load_si128(bam_aux + i);

		//Compare sequences
		misses += sse2_epi8_compare_and_count_diffs(v_ref, v_seq, &v_comp);

		//Store comparation values
		_mm_store_si128(comp_aux + i, v_comp);
	}

	//Set result
	memcpy(comp_res, comp_aux, bam_seq_l);

	//Set miss count
	if(miss_count)
	{
		*miss_count = misses;
	}

#ifdef _REENTRANT
	//Free
	_mm_free(ref_aux);
	_mm_free(bam_aux);
	_mm_free(comp_aux);
#endif

#else //Sequential Block

	//Iterates nucleotides in this read
	count = 0;
	for(i = 0; i < bam_seq_l; i++)
	{
		//0x00 Equals, 0x01 Diff
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
inline ERROR_CODE
nucleotide_miss_qual_sum(const char *ref_seq, const char *bam_seq, const char *bam_qual, size_t bam_seq_l, char *comp_seq, uint32_t *out_miss_count, uint32_t *out_sum_quals)
{
	uint32_t sum;
	int z;

	//Get missmatches
	nucleotide_compare(ref_seq, bam_seq, bam_seq_l, comp_seq, out_miss_count);

	//If wants quality sum
	if(out_sum_quals)
	{
		//Calculate miss sum
		sum = 0;
		if(*out_miss_count != 0)
		{
			for(z = 0; z < bam_seq_l; z++)
			{
				if(comp_seq[z] != 0)
					sum += bam_qual[z];
			}

			//Set output
			*out_sum_quals = sum;
		}
		else
		{
			//Perfect match
			*out_sum_quals = 0;
		}
	}

	return NO_ERROR;
}
