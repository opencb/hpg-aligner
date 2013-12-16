#ifndef BS_WRITER_H
#define BS_WRITER_H

#include <stdio.h>

#include "commons/log.h"
#include "commons/file_utils.h"

#include "commons/workflow_scheduler.h"

#include "bioformats/fastq/fastq_batch_reader.h"
#include "bioformats/bam/bam_file.h"


#include "buffers.h"
#include "batch_writer.h"

int bs_writer(void *data);

#endif // BS_WRITER_H
