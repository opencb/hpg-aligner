#include "aux_nucleotide.h"

ERROR_CODE
nucleotide_compare(char *ref_seq, char *bam_seq, size_t bam_seq_l, char *comp_res, uint32_t *miss_count)
{
	int i;
	uint32_t count;

	assert(ref_seq);
	assert(bam_seq);
	assert(bam_seq_l > 0);
	assert(comp_res);
	assert(miss_count);

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
}
