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

#include "aux_bam.h"
#include <assert.h>

ERROR_CODE
compare_bams_qual(const char* bamPath0, const char* bamPath1, const int cycles)
{
	bam_file_t* bamFile0;
	bam_file_t* bamFile1;
	bam_batch_t* bamBatch0;
	bam_batch_t* bamBatch1;
	bam1_t* bamAlig;
	alignment_t* aligAlig0;
	alignment_t* aligAlig1;
	int diff, i;

	printf("Opening BAM 1 form \"%s\" ...\n", bamPath0);
	printf("Opening BAM 2 form \"%s\" ...\n", bamPath1);
	bamFile0 = bam_fopen(bamPath0);
	bamFile1 = bam_fopen(bamPath1);
	printf("BAM opened!...\n");


	printf("\n\n---------------------------------------------------------\n");

	bamBatch0 = bam_batch_new(1, MULTIPLE_CHROM_BATCH);
	bamBatch1 = bam_batch_new(1, MULTIPLE_CHROM_BATCH);
	bam_fread_max_size(bamBatch0, 1, 1, bamFile0);
	bam_fread_max_size(bamBatch1, 1, 1, bamFile1);

	//Obtain first alignment from first bam
	bamAlig = bamBatch0->alignments_p[0];
	aligAlig0 = alignment_new_by_bam(bamAlig, 1);
	alignment_print(aligAlig0);

	//Obtain first alignment from second bam
	bamAlig = bamBatch1->alignments_p[0];
	aligAlig1 = alignment_new_by_bam(bamAlig, 1);
	alignment_print(aligAlig1);

	//Obtain quality diffs
	printf("Diffs: \nNuc\tQ1\tQ2\n");
	diff=0;
	for(i=0; i < 76; i++)
	{
		printf("%c \t%d ", aligAlig0->sequence[i], aligAlig0->quality[i]);
		if(aligAlig0->quality[i] == aligAlig1->quality[i])
		{
			printf("====\t%d\n", aligAlig1->quality[i]);
		}
		else
		{
			printf("\t%d\n", aligAlig1->quality[i]);
		}


		diff += abs(aligAlig1->quality[i] - aligAlig0->quality[i]);
	}
	printf("Total diff: %d\n", diff);

	printf("\n---------------------------------------------------------\n");
	printf("Closing BAMs...\n");
	bam_fclose(bamFile0);
	bam_fclose(bamFile1);
	bam_batch_free(bamBatch0, 1);
	bam_batch_free(bamBatch1, 1);

	printf("BAM closed.\n");

	return NO_ERROR;
}

ERROR_CODE
init_empty_bam_header(const unsigned int num_chroms, bam_header_t *header)
{
	int i;

	if(num_chroms == 0)
		return INVALID_INPUT_PARAMS_0;

	if(!header)
	{
		return INVALID_INPUT_PARAMS_NULL;
	}

	//Create a header with chroms targets number
	header->n_targets = num_chroms;
	header->target_name = (char **) calloc(num_chroms, sizeof(char *));
	header->target_len = (uint32_t*) calloc(num_chroms, sizeof(uint32_t));

	for(i = 0; i < num_chroms; i++)
	{
		header->target_name[i] = strdup("chr1");
		header->target_len[i] = strlen("chr1")+1;
	}

	header->text = strdup("@PG\tID:HPG-RECALIBRATOR\tVN:0.1\n");
	header->l_text = strlen(header->text);

	return NO_ERROR;
}

/**
 * Get string containing bam1 sequence nucleotides.
 */
char *
new_sequence_from_bam(bam1_t *bam1)
{
	char *seq;
	char *bam_seq = (char *)bam1_seq(bam1);
	int seq_len = bam1->core.l_qseq;

#ifdef __SSE2__
	seq = (char *) _mm_malloc(seq_len * sizeof(char), MEM_ALIG_SSE_SIZE);
#else
	seq = (char *) malloc(seq_len * sizeof(char));
#endif

	// nucleotide content
	for (int i = 0; i < seq_len; i++) {
		switch (bam1_seqi(bam_seq, i))
		{
		case 1:
			seq[i] = 'A';
			break;
		case 2:
			seq[i] = 'C';
			break;
		case 4:
			seq[i] = 'G';
			break;
		case 8:
			seq[i] = 'T';
			break;
		case 15:
			seq[i] = 'N';
			//printf("N");
			break;
		default:
			seq[i] = 'N';
			break;
		}
	}

	return seq;
}

/**
 * Get string containing bam1 sequence nucleotides.
 */
ERROR_CODE
new_sequence_from_bam_ref(bam1_t *bam1, char *seq, uint32_t max_l)
{
	char *bam_seq = (char *)bam1_seq(bam1);
	int seq_len = bam1->core.l_qseq;
	int i;

	if(seq_len > max_l)
		seq_len = max_l;

	// nucleotide content
	for (i = 0; i < seq_len; i++) {
		switch (bam1_seqi(bam_seq, i))
		{
		case 1:
			seq[i] = 'A';
			break;
		case 2:
			seq[i] = 'C';
			break;
		case 4:
			seq[i] = 'G';
			break;
		case 8:
			seq[i] = 'T';
			break;
		case 15:
			seq[i] = 'N';
			//printf("N");
			break;
		default:
			seq[i] = 'N';
			break;
		}
	}

	if(max_l > seq_len)
		seq[i] = '\0';

	return NO_ERROR;
}

/**
 * Get string containing bam1 quality for every nucleotide.
 */
char *
new_quality_from_bam(bam1_t *bam1, int base_quality)
{
	char *qual;
	char *bam_qual = (char *)bam1_qual(bam1);
	int qual_len = bam1->core.l_qseq;

#ifdef __SSE2__
	qual = (char *) _mm_malloc(qual_len * sizeof(char), MEM_ALIG_SSE_SIZE);
#else
	qual = (char *) malloc(qual_len * sizeof(char));
#endif
	for (int i = 0; i < qual_len; i++) {
		qual[i] = base_quality + bam_qual[i];
	}

	return qual;
}

/**
 * Get string containing bam1 quality for every nucleotide.
 */
ERROR_CODE
new_quality_from_bam_ref(bam1_t *bam1, U_QUALS base_quality, char *qual, U_CYCLES max_l)
{
	char *bam_qual = bam1_qual(bam1);
	int qual_len = bam1->core.l_qseq;

	if(qual_len > max_l)
		qual_len = max_l;

	for (int i = 0; i < qual_len; i++) {
		qual[i] = base_quality + bam_qual[i];
	}

	return NO_ERROR;
}

/**
 * Decomponse a string cigar in two vectors, containing number of ocurrences for every cigar operation (I, D, =, M, etc..)
 */
ERROR_CODE
decompose_cigar(char *cigar, uint8_t cigar_l, char *n_elem, char *type, uint8_t *types_l, uint8_t max_types_length)
{
	//uint8_t u_elems;
	//char c_type;
	char cigar_elem;
	int pos;
	int i;

	if(cigar == NULL
			|| cigar_l <= 0
			|| n_elem == NULL
			|| type == NULL
			|| types_l == NULL)
	{
		return INVALID_INPUT_PARAMS_NULL;
	}

	pos = 0;
	n_elem[pos] = 0;
	for(i = 0; i < cigar_l; i++)
	{
		cigar_elem = cigar[i];

		//Get number of elem
		if(cigar_elem <= '9' && cigar_elem >= '0')
		{
			//Is number
			n_elem[pos] *= 10;
			n_elem[pos] += cigar_elem - '0';
		}
		else
		{
			//Is type
			//switch(cigar_elem)
			//{
			//case 'I':
			//case 'D':
			//default:	//Missmatch, etc...
				if(n_elem[pos] > 0)	//Avoid void types like '0D', '0M' ..
				{
					type[pos] = cigar_elem;
					pos++;

					if(pos >= max_types_length)
					{
						//Reached max number of cigar components
						goto decompose_cigar_end;
					}

					n_elem[pos] = 0;
				}
				//break;
			//}
		}
	}

	//Set type length
	decompose_cigar_end:
	if(types_l != NULL)
		*types_l = pos;

	return NO_ERROR;
}

/**
 * Generates sequence with indels supressed. Insertions are deleted, Deletions are substituted by 'X'.
 */
ERROR_CODE
supress_indels(char *seq, U_CYCLES seq_l, char *cigar_elem, char *cigar_type, uint8_t cigar_type_l, char *seq_res, U_CYCLES *seq_res_l)
{
	int i, j, seq_i, res_i;
	char count;

	//Check null parameters
	if(seq == NULL
		|| cigar_elem == NULL
		|| cigar_type == NULL
		|| cigar_type_l == NULL
		|| seq_res == NULL
		|| seq_res_l == NULL)
	{
		return INVALID_INPUT_PARAMS_NULL;
	}

	//Iterate cigar
	for(i = 0, res_i = 0, seq_i = 0; i < cigar_type_l; i++)
	{
		switch(cigar_type[i])
		{
		case 'I':	//Insertion
			seq_i += cigar_elem[i];
			break;
		case 'D':	//Deletion
			count = cigar_elem[i];
			for(j = 0; j < count; j++)
			{
				seq_res[res_i] = 'X';
				res_i++;
			}
			break;
		default:	//Missmatch, etc...
			count = cigar_elem[i];
			memcpy(&seq_res[res_i], &seq[seq_i], count);
			res_i += count;
			seq_i += count;
			break;
		}
	}

	//Ser result length
	*seq_res_l = res_i;

	//Set string null character in las position
	seq_res[res_i] = '\0';

	return NO_ERROR;
}

/**
 * Get sequence and quality vector supressing insertions and adding deletions.
 */
ERROR_CODE
supress_indels_from_32_cigar(char *seq, char *qual, int32_t seq_l, uint32_t *cigar, uint16_t cigar_l, char *seq_res, char *qual_res, uint32_t *seq_res_l, uint32_t max_res_l)
{
	int i, j, seq_i, res_i;
	uint32_t count;
	uint32_t elem;
	uint8_t type;

	//Check null parameters
	if(seq == NULL
		|| cigar == NULL
		|| seq_res == NULL
		|| seq_res_l == NULL)
	{
		return INVALID_INPUT_PARAMS_NULL;
	}

	seq_i = 0;
	res_i = 0;
	for(i = 0; i < cigar_l; i++)
	{
		elem = cigar[i];
		count = elem >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		type = elem & BAM_CIGAR_MASK;	//Get type from cigar

		switch(type)
		{
		case BAM_CSOFT_CLIP:
		case BAM_CINS:	//Insertion
			if(count + seq_i > seq_l)
				count = seq_l - seq_i;
			seq_i += count;
			break;
		case BAM_CDEL:	//Deletion
			for(j = 0; j < count; j++)
			{
				seq_res[res_i] = 'X';
				if(qual && qual_res) qual_res[res_i] = '!';
				res_i++;
				if(res_i >= max_res_l - 1)	//Check array limit
				{
					res_i = max_res_l - 1;
					i = cigar_l;	//End loop
				}
			}
			break;
		case BAM_CMATCH:
		case BAM_CEQUAL:
		case BAM_CDIFF:
		case BAM_CHARD_CLIP:
		case BAM_CPAD:
			if(count + seq_i > seq_l)
				count = seq_l - seq_i;

			if(res_i + count >= max_res_l - 1)	//Check array limit
			{
				count = max_res_l - 1 - res_i;
				i = cigar_l;	//End loop
			}

			memcpy(&seq_res[res_i], &seq[seq_i], count * sizeof(char));
			if(qual && qual_res) memcpy(&qual_res[res_i], &qual[seq_i], count * sizeof(char));
			res_i += count;
			seq_i += count;
			break;

		default:
			fprintf(stderr, "WARNING: Unrecognised cigar N:%d T:%d\n", count, type);
			fflush(stderr);
			//abort();
		}
	}

	//Ser result length
	*seq_res_l = res_i;

	//Set string null character in last position
	seq_res[res_i] = '\0';
	qual_res[res_i] = '\0';

	return NO_ERROR;
}

ERROR_CODE
batch_count_chroms(bam_batch_t *batch, size_t *chrom_l)
{
	int i;
	int32_t last_chrom;
	int32_t actual_chrom;
	size_t result;

	//CHECK PARAMS
	{
		if(!batch || !chrom_l)
			return INVALID_INPUT_PARAMS_NULL;
	}

	last_chrom = -1;
	result = 0;
	for(i = 0; i < batch->num_alignments; i++)
	{
		//Set actual chrom
		actual_chrom = batch->alignments_p[i]->core.tid;

		//Count chroms
		if(actual_chrom != last_chrom)
		{
			last_chrom = actual_chrom;
			result++;
		}
	}

	//Set result
	*chrom_l = result;

	return NO_ERROR;
}

ERROR_CODE
batch_split_by_chrom(bam_batch_t *batch, bam_batch_t *v_batchs, size_t *res_batch_l, size_t max_res_batchs)
{
	int32_t last_chrom;
	int32_t actual_chrom;
	int i;

	bam_batch_t *actual_batch;
	int actual_batch_i;
	size_t count;

	//CHECK PARAMS
	{
		if(!batch || !v_batchs || !res_batch_l || max_res_batchs == 0)
			return INVALID_INPUT_PARAMS_NULL;
	}

	//Fill batchs
	actual_batch_i = 0;
	actual_batch = v_batchs;
	actual_chrom = -1;
	last_chrom = batch->alignments_p[0]->core.tid;
	count = 0;
	for(i = 0; i < batch->num_alignments; i++)
	{
		//Get next alignment chrom id
		actual_chrom = batch->alignments_p[i]->core.tid;

		//Check alignment chrom
		if(actual_chrom != last_chrom)
		{
			//Check valid count
			assert(count != 0);

			//Copy alignments to new batch
			actual_batch->alignments_p = (bam1_t **) malloc(sizeof(bam1_t *) * count);
			actual_batch->num_alignments = count;
			actual_batch->allocated_alignments = count;
			memcpy(actual_batch->alignments_p, &batch->alignments_p[i - count], count * sizeof(bam1_t*));
			//memset(&batch->alignments_p[i - count - 1], NULL, sizeof(bam1_t *) * count);	/* Set NULL in original batch */
			count = 0;

			//Set last chrom
			last_chrom = actual_chrom;
			actual_batch_i++;
			assert(actual_batch_i < max_res_batchs);	/* This must not happen */

			//Initialize batch
			actual_batch = &v_batchs[actual_batch_i];
			actual_batch->type = SINGLE_CHROM_BATCH;
			actual_batch->allocated_alignments = 0;
			actual_batch->num_alignments = 0;
		}

		//Count alignments for this batch
		count++;
	}

	//Copy alignments to new batch
	actual_batch->alignments_p = (bam1_t **) malloc(sizeof(bam1_t *) * count);
	actual_batch->num_alignments = count;
	actual_batch->allocated_alignments = count;
	memcpy(actual_batch->alignments_p, &batch->alignments_p[i - count], count * sizeof(bam1_t*));
	//memset(&batch->alignments_p[i - count - 1], NULL, sizeof(bam1_t *) * count);	/* Set NULL in original batch */

	//Set number of result batchs
	*res_batch_l = actual_batch_i + 1;

	return NO_ERROR;
}
