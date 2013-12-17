#ifndef PAIR_SERVER_H
#define PAIR_SERVER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "commons/commons.h"
#include "commons/system_utils.h"
#include "containers/list.h"

#include "buffers.h"
#include "timing.h"
#include "sw_server.h"

//====================================================================================
//  structures and prototypes
//====================================================================================

struct pair_server_input {
  pair_mng_t *pair_mng;
  report_optarg_t *report_optarg;

  list_t* pair_list;
  list_t* sw_list;
  list_t* write_list;
};

//------------------------------------------------------------------------------------

inline void generate_alignment_len(alignment_t *alig) {
  alig->map_len = 0;
  char *cigar_str = alig->cigar;
  int cigar_len = strlen(cigar_str);	  
  int c = 0;
  char op;
  char op_value[1024]; 

  for (int j = 0; j < cigar_len; j++) {
    op = cigar_str[j];
    if (op < 58) {
      op_value[c++] = op;
    } else {
      if (op == 'N' || op == 'M' || op == 'D') {
	op_value[c] = '\0';
	alig->map_len += atoi(op_value);
      }
      c = 0;
    }
  } 

}


void pair_server_input_init(pair_mng_t *pair_mng, report_optarg_t *report_optarg,
			    list_t* pair_list, list_t *sw_list,
			    list_t *write_list, pair_server_input_t* input);

//====================================================================================

void pair_server(pair_server_input_t* input);
void prepare_pair_server(pair_server_input_t* input);
//------------------------------------------------------------------------------------

int apply_pair(pair_server_input_t* input, batch_t *batch);

int prepare_alignments(pair_server_input_t* input, batch_t *batch);
int prepare_alignments_bs(pair_server_input_t* input, batch_t *batch);

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

#endif // PAIR_SERVER_H
