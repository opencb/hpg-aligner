#ifndef DATA_COLLECTOR_H_
#define DATA_COLLECTOR_H_

#include <stdlib.h>
#include <string.h>

#include <bioformats/bam/bam_file.h>
#include <aligners/bwt/genome.h>
#include "recalibrate/recal_config.h"
#include "aux/timestats.h"
#include "recalibrate/recal_config.h"
#include "bam_recal_library.h"

#include <omp.h>


#include "bioformats/bam/alignment.h"
#include "recal_structs.h"

long int unmapped;
long int duplicated;
long int notprimary;
#ifdef NOT_MAPPING_QUAL_ZERO
long int mapzero;
#endif
/***********************************************
 * BAM RECALIBRATION PHASE 1 - DATA COLLECT
 **********************************************/

/**
 * Private functions
 */

#endif /* DATA_COLLECTOR_H_ */
