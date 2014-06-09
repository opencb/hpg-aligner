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

#include "recal_structs.h"

#include <math.h>
#include <float.h>
#include <fenv.h>

/**
 *
 * DATA MANAGEMENT
 *
 */

/**
 * Initialize empty recalibration data struct.
 */
ERROR_CODE
recal_init_info(const U_CYCLES cycles, recal_info_t *out_data)
{
	recal_info_t *data;
	int vector_size;

	assert(out_data);

	vector_size = MAX_QUALITY - MIN_QUALITY;
	data = out_data;

	//Struct initialization
	data->min_qual = MIN_QUALITY;
	data->num_quals = vector_size;
	data->num_cycles = cycles;
	data->num_dinuc = NUM_DINUC;

	//Total counters
	data->total_miss = 0.0;
	data->total_bases = 0;
	data->total_delta = 0.0;
	data->total_estimated_Q = 0.0;

	//Quality vectors
	new_vector_double(vector_size, 0.0, &(data->qual_miss));
	data->qual_bases = (U_BASES *)malloc(vector_size * sizeof(U_BASES));
	new_vector_double(vector_size, 0.0, &(data->qual_delta));

	//Qual-Cycle matrix
	new_vector_double(vector_size * cycles, 0.0, &(data->qual_cycle_miss));
	data->qual_cycle_bases = (U_BASES *)malloc(vector_size * cycles * sizeof(U_BASES));
	new_vector_double(vector_size * cycles, 0.0, &(data->qual_cycle_delta));

	//Qual-Dinuc matrix
	new_vector_double(vector_size * NUM_DINUC, 0.0, &(data->qual_dinuc_miss));
	data->qual_dinuc_bases = (U_BASES *)malloc(vector_size * NUM_DINUC * sizeof(U_BASES));
	new_vector_double(vector_size * NUM_DINUC, 0.0, &(data->qual_dinuc_delta));

	//Initialize vector values
	memset(data->qual_bases, 0, vector_size * sizeof(U_BASES));
	memset(data->qual_cycle_bases, 0, vector_size * cycles * sizeof(U_BASES));
	memset(data->qual_dinuc_bases, 0, vector_size * NUM_DINUC * sizeof(U_BASES));

	return NO_ERROR;
}

/**
 * Free all resources of recalibration.
 */
ERROR_CODE
recal_destroy_info(recal_info_t *data)
{
	recal_info_t *d = data;

	//Free quality vector
	free(d->qual_miss);
	free(d->qual_bases);
	free(d->qual_delta);

	//Free quality-cycle matrix
	free(d->qual_cycle_miss);
	free(d->qual_cycle_bases);
	free(d->qual_cycle_delta);

	//Free quality-dinuc matrix
	free(d->qual_dinuc_miss);
	free(d->qual_dinuc_bases);
	free(d->qual_dinuc_delta);

	return NO_ERROR;
}

/**
 * Reduce data.
 */
ERROR_CODE
recal_reduce_info(recal_info_t *dst_data, recal_info_t *src_data)
{
	int i;

	//CHECK PARAMS
	{
		if(!dst_data || !src_data)
			return INVALID_INPUT_PARAMS_NULL;

		if(dst_data->num_cycles != src_data->num_cycles)
			return INVALID_INPUT_PARAMS;

		if(dst_data->min_qual != src_data->min_qual)
			return INVALID_INPUT_PARAMS;

		if(dst_data->num_quals != src_data->num_quals)
			return INVALID_INPUT_PARAMS;

		if(dst_data->num_dinuc != src_data->num_dinuc)
			return INVALID_INPUT_PARAMS;
	}

	//MERGE
	//TODO : Potential SSE implementation
	{
		//Total counters
		dst_data->total_miss += src_data->total_miss;
		dst_data->total_bases += src_data->total_bases;

		//Quality vectors
		for(i = 0; i < dst_data->num_quals; i++)
		{
			dst_data->qual_miss[i] += src_data->qual_miss[i];
			dst_data->qual_bases[i] += src_data->qual_bases[i];
		}

		//Qual-Cycle matrix
		for(i = 0; i < dst_data->num_quals * dst_data->num_cycles; i++)
		{
			dst_data->qual_cycle_miss[i] += src_data->qual_cycle_miss[i];
			dst_data->qual_cycle_bases[i] += src_data->qual_cycle_bases[i];
		}

		//Qual-Dinuc matrix
		for(i = 0; i < dst_data->num_quals * dst_data->num_dinuc; i++)
		{
			dst_data->qual_dinuc_miss[i] += src_data->qual_dinuc_miss[i];
			dst_data->qual_dinuc_bases[i] += src_data->qual_dinuc_bases[i];
		}
	}

	return NO_ERROR;
}

/**
 * Initialize empty data collection environment struct.
 */
ERROR_CODE
recal_get_data_init_env(const U_CYCLES cycles, recal_data_collect_env_t *collect_env)
{
	if(!collect_env)
		return INVALID_INPUT_PARAMS_NULL;

	//Sequence storage
	collect_env->bam_seq = (char *) _mm_malloc (sizeof(char) * cycles, MEM_ALIG_SSE_SIZE);
	collect_env->bam_quals = (char *) _mm_malloc (sizeof(char) * cycles, MEM_ALIG_SSE_SIZE);

	//Auxiliar
	collect_env->aux_res_seq = (char *) malloc(sizeof(char) * cycles);
	collect_env->aux_res_qual = (char *) malloc(sizeof(char) * cycles);

	//Maximum length
	collect_env->bam_seq_max_l = cycles;

	return NO_ERROR;
}

/**
 * Free all resources of data collect environment.
 */
ERROR_CODE
recal_get_data_destroy_env(recal_data_collect_env_t *collect_env)
{
	if(!collect_env)
		return INVALID_INPUT_PARAMS_NULL;

	_mm_free(collect_env->bam_seq);
	_mm_free(collect_env->bam_quals);

	free(collect_env->aux_res_seq);
	free(collect_env->aux_res_qual);

	free(collect_env);

	return NO_ERROR;
}

/**
 * Initialize empty recalibration environment struct.
 */
ERROR_CODE
recal_recalibration_init_env(const U_CYCLES cycles, recal_recalibration_env_t *recalibration_env)
{
	if(!recalibration_env)
		return INVALID_INPUT_PARAMS_NULL;

	//Quality storage
	recalibration_env->bam_quals = (char *) malloc (sizeof(char) * cycles);

	//Maximum length
	recalibration_env->bam_seq_max_l = cycles;

	return NO_ERROR;
}

/**
 * Free all resources of recalibration environment.
 */
ERROR_CODE
recal_recalibration_destroy_env(recal_recalibration_env_t *recalibration_env)
{
	if(!recalibration_env)
		return INVALID_INPUT_PARAMS_NULL;

	free(recalibration_env->bam_quals);


	free(recalibration_env);

	return NO_ERROR;
}

/**
 * Add recalibration data from one base.
 */
ERROR_CODE
recal_add_base(recal_info_t *data, const char qual, const U_CYCLES cycle, const char dinuc, const double match)
{
	int qual_index = qual - data->min_qual;
	int qual_cycle_index = qual_index * data->num_cycles + cycle;
	int qual_dinuc_index = qual_index * data->num_dinuc + dinuc;

	if(qual < MIN_QUALITY_TO_STAT)
		return NO_ERROR;

	//Error check
	if(qual_index >= MAX_QUALITY - MIN_QUALITY || qual_index < 0)
	{
		printf("add_base: ERROR, qual must be positive and minor than MAX_QUALITY - MIN_QUALITY ==> Qual = %d, Cycle = %d, Dinuc = %d, %s\n", qual, cycle, dinuc, match ? "Miss" : "");
		return INVALID_INPUT_QUAL;
	}
	if(cycle < 0 || cycle > data->num_cycles - 1)
	{
		printf("add_base: ERROR, cycle must be positive and minor than NUM_CYCLES (%d), Cycle: %d\n", data->num_cycles, cycle);
		return INVALID_INPUT_QUAL;
	}
	if(dinuc > NUM_DINUC - 1)
	{
		printf("add_base: ERROR, dinuc must be minor than NUM_DINUC %d \n", dinuc);
		return INVALID_INPUT_DINUC;
	}

	//Increase total counters
	if(!match)
	data->total_miss++;
	data->total_bases++;

	//Increase quality vector
	if(!match)
	data->qual_miss[qual_index]++;
	data->qual_bases[qual_index]++;

	//Increase quality-cycle matrix
	if(!match)
	data->qual_cycle_miss[qual_cycle_index]++;
	data->qual_cycle_bases[qual_cycle_index]++;

	//Increase quality-dinuc matrix (if not dinuc != -1)
	if(dinuc >= 0)
	{
		if(!match)
		data->qual_dinuc_miss[qual_dinuc_index]++;
		data->qual_dinuc_bases[qual_dinuc_index]++;
	}
	else
	{
		printf("ERROR: unrecognized dinuc Q: %d, C: %d, D: %d, Match: %d\n", qual, cycle, dinuc, match);
		return INVALID_INPUT_DINUC;
	}

	return NO_ERROR;
}

/**
 * Add recalibration data from vector of bases
 */
ERROR_CODE
recal_add_base_v(recal_info_t *data, const char *seq, const char *quals, const U_CYCLES init_cycle, const U_CYCLES num_cycles, const char *dinuc, const char *matches)
{
	U_CYCLES i;
	U_CYCLES cycles;
	ERROR_CODE err;

	cycles = num_cycles;

	//Iterates cycles
	//for(i = init_cycle; i <= end_cycle; i++)
	for(i = 0; i < cycles; i++)
	{
		switch(seq[i])
		{
		case 'A':
		case 'C':
		case 'G':
		case 'T':
		#ifndef NOT_COUNT_NUCLEOTIDE_N
		case 'N':
		#endif
			err = recal_add_base(data, quals[i], i + init_cycle, dinuc[i], (double)matches[i]);
			if(err)
			{
				switch(err)
				{
				case INVALID_INPUT_QUAL:
					printf("Error INVALID_INPUT_QUAL Base: %c\n", seq[i]);
					break;
				case INVALID_INPUT_DINUC:
					printf("Error INVALID_INPUT_DINUC Base: %c\n", seq[i]);
					break;
				default:
					printf("Error UNKNOWN Base: %c\n", seq[i]);
					break;
				}
			}
			break;

		default:
			//printf("ERROR: Corrupted nucleotide read = %c\n", seq[i]);
			break;
		}
	}

	return NO_ERROR;
}

/**
 * Compute deltas from bases and misses.
 */
ERROR_CODE
recal_calc_deltas(recal_info_t* data)
{
	double phred, err0;
	//double r_empirical;
	int matrix_index;
	int i, j;

	double estimated_Q;
	double est2;
	double emp_Q;
	double delta;

	//Time measures
	#ifdef D_TIME_DEBUG
		time_init_slot(D_SLOT_PROCCESS_DELTAS, TIME_GLOBAL_STATS);
	#endif

	printf("Processing deltas...\n");

	//Estimated Q
	recal_get_estimated_Q(data->qual_bases, data->num_quals, (U_QUALS)0, &estimated_Q);

	//Get empirical global qual
	recal_get_empirical_Q(data->total_miss, data->total_bases, estimated_Q, &emp_Q);

	//Calc global delta
	data->total_estimated_Q = estimated_Q;
	data->total_delta = emp_Q - estimated_Q;
	printf("Estimated %.2f \tEmpirical %.2f \t TotalDelta %.2f\n", estimated_Q, emp_Q, data->total_delta);

	//Delta R
	for(i = 0; i < data->num_quals; i++)
	{
		if(data->qual_bases[i] != 0)
		{
			//int calidad = Qvalue( (double)(data->qual_miss[i]) / (double)(data->qual_bases[i]) );
			//data->qual_delta[i] = calidad
			//		- (double)(i /*+ data->min_qual*/)
			//		- data->total_delta;

			recal_get_empirical_Q(data->qual_miss[i], data->qual_bases[i], /*i*/data->total_delta + estimated_Q, &emp_Q);
			delta = emp_Q - (data->total_delta + estimated_Q);
			data->qual_delta[i] = delta;
			//printf("Qual %d  \tEmpirical %.2f  \tDelta %.2f\n", i, emp_Q, data->qual_delta[i]);
		}
		else
		{
			data->qual_delta[i] = 0.0;
		}
	}


	//Delta R,C
	for(i = 0; i < data->num_quals; i++)
	{
		for(j = 0; j < data->num_cycles; j++)
		{
			matrix_index = i * data->num_cycles + j;
			if(data->qual_cycle_bases[matrix_index] != 0)
			{
				//data->qual_cycle_delta[matrix_index] = Qvalue((double)(data->qual_cycle_miss[matrix_index]) / (double)(data->qual_cycle_bases[matrix_index]))
				//	- (double)(i /*+ data->min_qual*/)
				//	- (data->total_delta + data->qual_delta[i]);
				recal_get_empirical_Q(data->qual_cycle_miss[matrix_index], data->qual_cycle_bases[matrix_index],
						data->qual_delta[i] + data->total_delta + estimated_Q, &emp_Q);
				delta = emp_Q - (data->qual_delta[i] + data->total_delta + estimated_Q);
				data->qual_cycle_delta[matrix_index] = delta;
			}
			else
			{
				data->qual_cycle_delta[matrix_index] = 0.0;
			}
		}
	}

	//Delta R,D
	for(i = 0; i < data->num_quals; i++)
	{
		for(j = 0; j < data->num_dinuc; j++)
		{
			matrix_index = i * data->num_dinuc + j;
			if(data->qual_dinuc_bases[matrix_index] != 0)
			{
				//data->qual_dinuc_delta[matrix_index] = Qvalue((double)(data->qual_dinuc_miss[matrix_index]) / (double)(data->qual_dinuc_bases[matrix_index]))
				//	- (double)(i /*+ data->min_qual*/)
				//	- (data->total_delta + data->qual_delta[i]);
				recal_get_empirical_Q(data->qual_dinuc_miss[matrix_index], data->qual_dinuc_bases[matrix_index],
						data->qual_delta[i] + data->total_delta + estimated_Q, &emp_Q);
				delta = emp_Q - (data->qual_delta[i] + data->total_delta + estimated_Q);
				data->qual_dinuc_delta[matrix_index] = delta;
			}
			else
			{
				data->qual_dinuc_delta[matrix_index] = 0.0;
			}
		}
	}

	printf("Deltas processed.\n");

	#ifdef D_TIME_DEBUG
		time_set_slot(D_SLOT_PROCCESS_DELTAS, TIME_GLOBAL_STATS);
	#endif

	return NO_ERROR;
}

/**
 *
 * ENUMERATION FUNCTIONS
 *
 */

/**
 * Return dinucleotide enumeration from two bases.
 */
ERROR_CODE
recal_get_dinuc(const char A, const char B, U_DINUC *out_dinuc)
{
	*out_dinuc = d_X;

	switch(A)
	{
			case 'A':
			switch(B)
			{
					case 'A':
						*out_dinuc = dAA;
					break;
					case 'G':
						*out_dinuc = dAG;
					break;
					case 'C':
						*out_dinuc = dAC;
					break;
					case 'T':
						*out_dinuc = dAT;
					break;
			}
			break;
			case 'G':
			switch(B)
			{
					case 'A':
						*out_dinuc = dGA;
					break;
					case 'G':
						*out_dinuc = dGG;
					break;
					case 'C':
						*out_dinuc = dGC;
					break;
					case 'T':
						*out_dinuc = dGT;
					break;
			}
			break;
			case 'C':
			switch(B)
			{
					case 'A':
						*out_dinuc = dCA;
					break;
					case 'G':
						*out_dinuc = dCG;
					break;
					case 'C':
						*out_dinuc = dCC;
					break;
					case 'T':
						*out_dinuc = dCT;
					break;
			}
			break;
			case 'T':
			switch(B)
			{
					case 'A':
						*out_dinuc = dTA;
					break;
					case 'G':
						*out_dinuc = dTG;
					break;
					case 'C':
						*out_dinuc = dTC;
					break;
					case 'T':
						*out_dinuc = dTT;
					break;
			}
			break;
			case '_':
			case 'N':

			break;
	}

	return NO_ERROR;
}

/**
 *
 * FILE OPERATIONS
 *
 */

/**
 * Print to file data from recalibration.
 */
ERROR_CODE
recal_fprint_info(const recal_info_t *data, const char *path)
{
	FILE *fp;
	int i,j;
	int n_quals = data->num_quals;
	int n_cycles = data->num_cycles;
	int n_dinuc = data->num_dinuc;

	fp = fopen(path, "w+");

	//Info msg
	printf("\n----------------\nPrinting data on \"%s\" file...\n----------------\n", path);

	//Print general info
	fprintf(fp, "==============================\nGENERAL: \nTotal miss: %.2lf\nTotal bases: %lu\nTotal delta: %.2lf\nEstimated Q: %.2lf\n",
			(double)data->total_miss, (U_BASES)data->total_bases, (double)data->total_delta, (double)data->total_estimated_Q);

	//Print quality infos
	fprintf(fp, "==============================\nQUAL VECTOR:\n");
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%3d %8.2f %10u %6.2f \n", i + MIN_QUALITY, data->qual_miss[i] /*- SMOOTH_CONSTANT_MISS*/, data->qual_bases[i] /*- SMOOTH_CONSTANT_BASES*/, data->qual_delta[i] /*- SMOOTH_CONSTANT_MISS*/);
	}

	//Print cycle infos
	fprintf(fp, "==============================\nQUAL-CYCLE MISS MATRIX (%d Quals - %d Cycles):\n", data->num_quals, data->num_cycles);
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%d \t", i + MIN_QUALITY);
		for(j = 0; j < n_cycles; j++)
		{
			fprintf(fp, "%.2f \t", data->qual_cycle_miss[i * n_cycles + j] /*- SMOOTH_CONSTANT_MISS*/);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "==============================\nQUAL-CYCLES BASES MATRIX (%d Quals - %d Cycles):\n", data->num_quals, data->num_cycles);
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%d \t", i + MIN_QUALITY);
		for(j = 0; j < n_cycles; j++)
		{
			fprintf(fp, "%u \t", data->qual_cycle_bases[i * n_cycles + j] /*- SMOOTH_CONSTANT_BASES*/);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "==============================\nQUAL-CYCLE DELTA MATRIX (%d Quals - %d Cycles):\n", data->num_quals, data->num_cycles);
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%d \t", i + MIN_QUALITY);
		for(j = 0; j < n_cycles; j++)
		{
			fprintf(fp, "%.2f \t", data->qual_cycle_delta[i * n_cycles + j]);
		}
		fprintf(fp, "\n");
	}

	//Print dinuc infos
	fprintf(fp, "==============================\nQUAL-DINUC MISS MATRIX (%d Quals - %d Dinuc):\n", data->num_quals, data->num_dinuc);
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%d \t", i + MIN_QUALITY);
		for(j = 0; j < n_dinuc; j++)
		{
			fprintf(fp, "%.2f \t", data->qual_dinuc_miss[i * n_dinuc + j] /*- SMOOTH_CONSTANT_MISS*/);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "==============================\nQUAL-DINUC BASES MATRIX (%d Quals - %d Dinuc):\n", data->num_quals, data->num_dinuc);
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%d \t", i + MIN_QUALITY);
		for(j = 0; j < n_dinuc; j++)
		{
			fprintf(fp, "%u \t", data->qual_dinuc_bases[i * n_dinuc + j] /*- SMOOTH_CONSTANT_BASES*/);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "==============================\nQUAL-DINUC DELTA MATRIX (%d Quals - %d Dinuc):\n", data->num_quals, data->num_dinuc);
	for(i = 0; i < n_quals; i++)
	{
		fprintf(fp, "%d \t", i + MIN_QUALITY);
		for(j = 0; j < n_dinuc; j++)
		{
			fprintf(fp, "%.2f \t", data->qual_dinuc_delta[i * n_dinuc + j]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);

	return NO_ERROR;
}

/**
 * Save to file recalibration data.
 */
ERROR_CODE
recal_save_recal_info(const recal_info_t *data, const char *path)
{
	FILE *fp;

	printf("\n----------------\nSaving recalibration data to \"%s\"\n----------------\n", path);

	fp = fopen(path, "w+");

	fwrite(data, sizeof(U_QUALS), 1, fp);	//min_qual
	fwrite(data, sizeof(U_QUALS), 1, fp);	//num_quals
	fwrite(data, sizeof(U_CYCLES), 1, fp);	//num_cycles
	fwrite(data, sizeof(U_DINUC), 1, fp);	//num_dinuc

	//Save total counters
	fwrite(&data->total_miss, sizeof(double), 1, fp);
	fwrite(&data->total_bases, sizeof(U_BASES), 1, fp);
	fwrite(&data->total_delta, sizeof(double), 1, fp);
	fwrite(&data->total_estimated_Q, sizeof(double), 1, fp);

	//Save qual counters
	fwrite(data->qual_miss, sizeof(double), data->num_quals, fp);
	fwrite(data->qual_bases, sizeof(U_BASES), data->num_quals, fp);
	fwrite(data->qual_delta, sizeof(double), data->num_quals, fp);

	//Save cycle counters
	fwrite(data->qual_cycle_miss, sizeof(double), data->num_quals * data->num_cycles, fp);
	fwrite(data->qual_cycle_bases, sizeof(U_BASES), data->num_quals * data->num_cycles, fp);
	fwrite(data->qual_cycle_delta, sizeof(double), data->num_quals * data->num_cycles, fp);

	//Save dinuc counters
	fwrite(data->qual_dinuc_miss, sizeof(double), data->num_quals * data->num_dinuc, fp);
	fwrite(data->qual_dinuc_bases, sizeof(U_BASES), data->num_quals * data->num_dinuc, fp);
	fwrite(data->qual_dinuc_delta, sizeof(double), data->num_quals * data->num_dinuc, fp);

	fclose(fp);

	return NO_ERROR;
}

/**
 * Load from file recalibration data.
 */
ERROR_CODE
recal_load_recal_info(const char *path, recal_info_t *data)
{
	FILE *fp;

	printf("\n----------------\nLoading recalibration data \"%s\"\n----------------\n", path);

	fp = fopen(path, "r");

	fread(data, sizeof(U_QUALS), 1, fp);	//min_qual
	fread(data, sizeof(U_QUALS), 1, fp);	//num_quals
	fread(data, sizeof(U_CYCLES), 1, fp);	//num_cycles
	fread(data, sizeof(U_DINUC), 1, fp);	//num_dinuc

	//Read total counters
	fread(&data->total_miss, sizeof(double), 1, fp);
	fread(&data->total_bases, sizeof(U_BASES), 1, fp);
	fread(&data->total_delta, sizeof(double), 1, fp);
	fread(&data->total_estimated_Q, sizeof(double), 1, fp);

	//Read qual counters
	fread(data->qual_miss, sizeof(double), data->num_quals, fp);
	fread(data->qual_bases, sizeof(U_BASES), data->num_quals, fp);
	fread(data->qual_delta, sizeof(double), data->num_quals, fp);

	//Read cycle counters
	fread(data->qual_cycle_miss, sizeof(double), data->num_quals * data->num_cycles, fp);
	fread(data->qual_cycle_bases, sizeof(U_BASES), data->num_quals * data->num_cycles, fp);
	fread(data->qual_cycle_delta, sizeof(double), data->num_quals * data->num_cycles, fp);

	//Read dinuc counters
	fread(data->qual_dinuc_miss, sizeof(double), data->num_quals * data->num_dinuc, fp);
	fread(data->qual_dinuc_bases, sizeof(U_BASES), data->num_quals * data->num_dinuc, fp);
	fread(data->qual_dinuc_delta, sizeof(double), data->num_quals * data->num_dinuc, fp);

	fclose(fp);

	return NO_ERROR;
}

/**
 * PRIVATE FUNCTIONS
 */
