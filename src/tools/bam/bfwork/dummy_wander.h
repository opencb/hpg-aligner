/*
 * dummy_wander.h
 *
 *  Created on: 20/06/2014
 *      Author: rmoreno
 */

#ifndef DUMMY_WANDER_H_
#define DUMMY_WANDER_H_

#include <assert.h>

#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <stddef.h>
#include <unistd.h>

#include <omp.h>

#include "aux/aux_common.h"
#include "aux/aux_library.h"
#include "aux/timestats.h"
#include "bfwork/bfwork.h"

EXTERNC ERROR_CODE dummy_bam_file(const char *bam_path, const char *ref_path, const char *outbam, const char *stats_path);

#endif /* DUMMY_WANDER_H_ */
