#ifndef FILTER_BAM_H
#define FILTER_BAM_H

/*
 * filter_bam.h
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include <stdlib.h>
#include <string.h>

//#include "argtable2.h"
//#include "libconfig.h"
//#include "commons/log.h"
//#include "commons/system_utils.h"
//#include "commons/file_utils.h"

#include "commons/workflow_scheduler.h"
#include "containers/array_list.h"
#include "containers/khash.h"

#include "bioformats/bam/bam_file.h"
#include "bioformats/bam/bam_filter.h"

#include "commons_bam.h"
#include "filter_options.h"

//------------------------------------------------------------------------

void filter_bam(filter_options_t *opts);

//------------------------------------------------------------------------

#endif // end of FILTER_BAM_H

//------------------------------------------------------------------------
//------------------------------------------------------------------------
