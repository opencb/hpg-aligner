#ifndef RECAL_STRUCTS_H_
#define RECAL_STRUCTS_H_

#include "aux/aux_common.h"
#include "aux/aux_library.h"
#include "aux/timestats.h"
#include "recalibrate/recal_config.h"
#include "bam_recal_library.h"


struct data_collect_env {

	//Sequence storage
	char *bam_seq;
	char *bam_quals;

	//Auxiliar
	char *aux_res_seq;
	char *aux_res_qual;

	//Maximum length
	U_CYCLES bam_seq_max_l;

};

struct recalibration_env {

	//Quality storage
	char *bam_quals;

	//Maximum length
	U_CYCLES bam_seq_max_l;

};


/**
 * PRIVATE FUNCTIONS
 */

#endif /* RECAL_STRUCTS_H_ */
