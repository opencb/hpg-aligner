#ifndef MULTI_ALIGNER_H
#define MULTI_ALIGNER_H

#include <mpi.h>

#include "buffers.h"
#include "options.h"
#include "breakpoint.h"
#include "rna/rna_splice.h"
#include "mpi/rna_aligner_mpi.h"

//For SA Mapping
#include "sa/sa_index3.h"
#include "sa/sa_search.h"

#include "rna/workflow_scheduler_SA.h"
#include "commons/workflow_scheduler.h"
#include "aligners/bwt/bwt.h"
#include "bioformats/bam/bam_file.h"
#include "bioformats/fastq/fastq_file.h"
#include "containers/array_list.h"
#include "containers/linked_list.h"

int hpg_multialigner_main(int argc, char *argv[]);



















#endif
