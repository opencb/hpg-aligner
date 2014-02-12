#include "aux_nucleotide.h"

/**
 * Compare two sequences and obtain missmatches.
 */
ERROR_CODE
nucleotide_compare(char *ref_seq, char *bam_seq, size_t bam_seq_l, char *comp_res, uint32_t *miss_count)
{
	int i;
	uint32_t count;

	assert(ref_seq);
	assert(bam_seq);
	assert(bam_seq_l > 0);
	assert(comp_res);
	//assert(miss_count);

	//SSE
	#ifdef __SSE2__
	__m128i v_ref, v_seq, v_comp;
	#endif

	//Iterates nucleotides in this read
	for(i = 0; i < bam_seq_l; i++)
	{
/*#ifdef __SSE2__	//SSE Block
		if( (i + 16) < bam_seq_l)
		{
			//Use SSE
			//_mm_prefetch(&ref_seq[i + 16], _MM_HINT_T0);
			//_mm_prefetch(&bam_seq[i + 16], _MM_HINT_T0);

			//Pack sequences
			v_ref = _mm_load_si128(&ref_seq[i]);
			v_seq = _mm_load_si128(&bam_seq[i]);

			//Compare sequences
			v_comp = _mm_cmpeq_epi8(v_ref, v_seq);

			//Store comparation values
			_mm_store_si128(&comp_res[i], v_comp);

			i += 15;
		}
		else
#endif //SSE Block*/
		{
			if(ref_seq[i] != bam_seq[i])
			{
				comp_res[i] = 0x00;	//Diff
			}
			else
			{
				comp_res[i] = 0xFF;	//Equals
			}
		}
	}

	if(miss_count)
	{
		count = 0;
		for(i = 0; i < bam_seq_l; i++)
		{
			if(!comp_res[i])
			{
				count++;
			}
		}
		*miss_count = count;
	}

	return NO_ERROR;
}

/**
 * Compare two sequences and obtain missmatches and missmatches qualities summatory.
 */
ERROR_CODE
nucleotide_miss_qual_sum(char *ref_seq, char *bam_seq, char *bam_qual, size_t bam_seq_l, char *comp_seq, uint32_t *out_miss_count, uint32_t *out_sum_quals)
{
	uint32_t sum;

	//Get missmatches
	nucleotide_compare(ref_seq, bam_seq, bam_seq_l, comp_seq, out_miss_count);

	//If wants quality sum
	if(out_sum_quals)
	{
		//Calculate miss sum
		sum = 0;
		if(*out_miss_count != 0)
		{
			int z;
			for(z = 0; z < bam_seq_l; z++)
			{
				if(comp_seq[z] == 0)
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
