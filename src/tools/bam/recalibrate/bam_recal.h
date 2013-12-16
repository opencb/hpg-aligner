#ifndef BAM_RECAL_H_
#define BAM_RECAL_H_

//#include <pthread.h>
#include <omp.h>

#include <bioformats/bam/bam_file.h>

#include "aux/aux_common.h"
#include "aux/aux_library.h"
#include "aux/timestats.h"
#include "recalibrate/recal_config.h"
#include "bam_recal_library.h"



#include "recal_structs.h"

/***********************************************
 * BAM RECALIBRATION PHASE 2 - RECALIBRATION
 **********************************************/

typedef struct batch_out{
	bam_batch_t *batch;
	recal_info_t *info;
} batch_out_t;

/**
 * Private functions
 */

static INLINE void recal_recalibrate_alignment_priv(const bam1_t* alig, const recal_info_t *bam_info, recal_recalibration_env_t *recalibration_env);


#endif /* BAM_RECAL_H_ */
