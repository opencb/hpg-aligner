#ifndef STATS_BAM_H
#define STATS_BAM_H

/*
 * stats_bam.h
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

#include "sqlite/sqlite3.h"
#include "samtools/bam.h"

#include "commons/workflow_scheduler.h"
#include "containers/khash.h"
#include "containers/array_list.h"
#include "bioformats/features/region/region_table.h"
#include "bioformats/db/db_utils.h"
#include "bioformats/bam/bam_file.h"
#include "bioformats/bam/bam_db.h"
#include "bioformats/bam/bam_stats.h"
#include "bioformats/bam/bam_filter.h"


#include "stats_options.h"
#include "stats_report.h"

//------------------------------------------------------------------------

void stats_bam(stats_options_t *opts);

//------------------------------------------------------------------------

#endif // end of STATS_BAM_H

//------------------------------------------------------------------------
//------------------------------------------------------------------------
