#ifndef DNA_ALIGNER_H
#define DNA_ALIGNER_H

#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

#include "commons/log.h"
#include "commons/file_utils.h"

#include "commons/workflow_scheduler.h"

#include "bioformats/fastq/fastq_batch_reader.h"

#include "error.h"
#include "timing.h"
#include "buffers.h"
#include "bwt_server.h"
#include "batch_writer.h"
#include "region_seeker.h"
#include "cal_seeker.h"
#include "sw_server.h"
#include "pair_server.h"
#include "options.h"
#include "statistics.h"
#include "workflow_functions.h"

#include "methylation.h"
#include "bs_writer.h"


//void run_bs_aligner(genome_t *genome, genome_t *genome1, genome_t *genome2,
void run_bs_aligner(genome_t *genome2, genome_t *genome1, genome_t *genome,
		    bwt_index_t *bwt_index2, bwt_index_t *bwt_index1, 
		    bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		    pair_mng_t *pair_mng, report_optarg_t *report_optarg,
		    options_t *options);

#endif
