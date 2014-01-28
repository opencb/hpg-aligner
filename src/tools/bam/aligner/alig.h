/*
 * alig.h
 *
 *  Created on: Dec 12, 2013
 *      Author: rmoreno
 */

#ifndef ALIG_H_
#define ALIG_H_

#include <assert.h>

#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <stddef.h>

#include <omp.h>

#include <bioformats/bam/samtools/bam.h>
#include <bioformats/bam/bam_file.h>
#include <aligners/bwt/genome.h>
#include "aux/aux_common.h"
#include "aux/aux_library.h"
#include "aux/timestats.h"
#include "alig_region.h"

//OPTIONS
#define ALIG_LEFT_ALIGN 0x01

#define ALIG_LIST_IN_SIZE	10000
#define ALIG_LIST_COUNT_THRESHOLD_TO_WRITE 1000

#define ALIG_REFERENCE_ADDITIONAL_OFFSET 100

/**
 * Time measures
 */
#define D_TIME_DEBUG
#ifdef D_TIME_DEBUG
	enum alig_slots {
		//GENERAL
		D_SLOT_TOTAL,
		D_SLOT_PROCCESS,
		D_SLOT_INIT,

		//BAM I/0
		D_SLOT_READ,
		D_SLOT_WRITE,

		//REALIGN
		D_SLOT_NEXT,
		D_SLOT_REFERENCE_LOAD,
		D_SLOT_HAPLO_GET,
		D_SLOT_REALIG_PER_HAPLO
	};

#endif


/**
 * REALIGNMENT CONTEXT
 */

typedef struct {
	char *reference;
	size_t position;
	size_t length;
} alig_reference_t;

typedef struct {
	//Input list
	linked_list_t *in_list;

	//Reference genome
	genome_t *genome;

	//BAM lists
	array_list_t *process_list;

	//Alignments readed
	size_t read_count;
	size_t last_readed_count;

	//Haplotypes
	array_list_t *haplo_list;

	//Current region
	alig_region_t region;

	//Current reference
	alig_reference_t reference;

	//Left align
	uint8_t flags;

} alig_context_t;


/**
 * INTERVAL STATUS
 */

typedef enum {
	NO_INTERVAL,
	INTERVAL
} alig_status;

/**
 * CONTEXT
 */

EXTERNC ERROR_CODE alig_init(alig_context_t *context, linked_list_t *in_list, genome_t *genome);
EXTERNC ERROR_CODE alig_destroy(alig_context_t *context);
EXTERNC ERROR_CODE alig_validate(alig_context_t *context);

/**
 * REGION OPERATIONS
 */

EXTERNC ERROR_CODE alig_region_next(alig_context_t *context);
EXTERNC ERROR_CODE alig_region_load_reference(alig_context_t *context);
EXTERNC ERROR_CODE alig_region_haplotype_process(alig_context_t *context);
EXTERNC ERROR_CODE alig_region_indel_realignment(alig_context_t *context);
EXTERNC ERROR_CODE alig_region_clear(alig_context_t *context);

EXTERNC ERROR_CODE alig_bam_file2(char *bam_path, char *ref_name, char *ref_path);

/**
 * BAM REALIGN
 */

EXTERNC ERROR_CODE alig_bam_list_realign(array_list_t *bam_list, array_list_t *haplotype_list, genome_t* ref);



#endif /* ALIG_H_ */
