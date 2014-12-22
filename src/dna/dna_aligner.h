#ifndef DNA_ALIGNER_H
#define DNA_ALIGNER_H

#include <stdio.h>
#include <stdlib.h>

#include "htslib/hts.h"

#include "commons/workflow_scheduler.h"
#include "bioformats/bam/bam_file.h"

#include "options.h"
#include "batch_writer.h"

#include "tools/bam/bfwork/bam_file_ops.h"

#include "sa/sa_index3.h"
#include "sa/sa_search.h"

#include "dna/sa_dna_commons.h"
#include "dna/doscadfun.h"
#include "dna/sa_io_stages.h"
#include "dna/sa_mapper_stage.h"


//--------------------------------------------------------------------

void dna_aligner(options_t *options);

//--------------------------------------------------------------------
//--------------------------------------------------------------------
#endif // DNA_ALIGNER_H
