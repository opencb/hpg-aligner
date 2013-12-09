#ifndef RECAL_STRUCTS_H_
#define RECAL_STRUCTS_H_

#include <recal_common.h>
#include <recal_config.h>
#include <aux_library.h>
#include <bam_recal_library.h>
#include <timestats.h>

/**
 * Recalibration data storage
 * This struct hold all data necessary to perform recalibrations
 */
struct recal_info {
	U_QUALS min_qual;	//Set minor stored quality
	U_QUALS num_quals;	//Set range of stored qualities
	U_CYCLES num_cycles;	//Set maximum number of cycles stored
	U_DINUC num_dinuc;	//Set number of maximum dinucleotides

	double total_miss;				//Total misses
	U_BASES total_bases;				//Total bases
	double total_delta;			//Global delta
	double total_estimated_Q;	//Global estimated Quality

	double* qual_miss;				//Misses per quality
	U_BASES* qual_bases;				//Bases per quality
    double* qual_delta;			//Delta per quality

    double* qual_cycle_miss;		//Misses per quality-cycle pair
    U_BASES* qual_cycle_bases;		//Bases per quality-cycle pair
    double* qual_cycle_delta;		//Deltas per quality-cycle pair

    double* qual_dinuc_miss;		//Misses per quality-dinuc pair
    U_BASES* qual_dinuc_bases;		//Bases per quality-dinuc pair
    double* qual_dinuc_delta;		//Deltas per quality-dinuc pair
};


typedef struct data_collect_env {

	//Sequence storage
	char *bam_seq;
	char *bam_quals;

	//Auxiliar
	char *aux_res_seq;
	char *aux_res_qual;

	//Maximum length
	U_CYCLES bam_seq_max_l;

};

typedef struct recalibration_env {

	//Quality storage
	char *bam_quals;

	//Maximum length
	U_CYCLES bam_seq_max_l;

};

/**
 * PRIVATE FUNCTIONS
 */

#endif /* RECAL_STRUCTS_H_ */
