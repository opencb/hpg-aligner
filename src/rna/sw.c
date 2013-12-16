#include "sw.h"

extern struct timeval t1, t2;
extern double t_diff, t_total; 

//====================================================================================
// Smith-Waterman structures and functions (SIMD version)
//====================================================================================

sw_simd_input_t* sw_simd_input_new(unsigned int depth) {
     sw_simd_input_t* input_p = calloc(1, sizeof(sw_simd_input_t));
     
     input_p->depth = depth;
     input_p->seq_p = (char**) calloc(depth, sizeof(char*));
     input_p->ref_p = (char**) calloc(depth, sizeof(char*));
     input_p->seq_len_p = (unsigned int*) calloc(depth, sizeof(unsigned int));
     input_p->ref_len_p = (unsigned int*) calloc(depth, sizeof(unsigned int));
     
     return input_p;
}

//------------------------------------------------------------------------------------

void sw_simd_input_free(sw_simd_input_t* input_p) {
     if (input_p == NULL) return;

     if (input_p->seq_p != NULL) free(input_p->seq_p);
     if (input_p->ref_p != NULL) free(input_p->ref_p);
     if (input_p->seq_len_p != NULL) free(input_p->seq_len_p);
     if (input_p->ref_len_p != NULL) free(input_p->ref_len_p);
     
     free(input_p);
}

//------------------------------------------------------------------------------------

void sw_simd_input_add(char* seq_p, unsigned int seq_len, char* ref_p,
		       unsigned int ref_len, unsigned int index, sw_simd_input_t* input_p) {
     if (index >= input_p->depth) {
	  printf("ERROR: out of range!\n");
	  return;
     }
     
     input_p->seq_p[index] = seq_p;
     input_p->ref_p[index] = ref_p;
     input_p->seq_len_p[index] = seq_len;
     input_p->ref_len_p[index] = ref_len;
}

//------------------------------------------------------------------------------------

void sw_simd_input_display(unsigned int depth, sw_simd_input_t* input_p) {
     char str[2048];
     
     printf("sw_simd_input (depth = %i of %i):\n", depth, input_p->depth);
     for(int i = 0; i < depth; i++) {	  
	  memcpy(str, input_p->seq_p[i], input_p->seq_len_p[i]);
	  str[input_p->seq_len_p[i]] = '\0';
	  printf("%i:\n", i);
	  printf("seq: %s (%i)\n", str, input_p->seq_len_p[i]);
	  
	  memcpy(str, input_p->ref_p[i], input_p->ref_len_p[i]);
	  str[input_p->ref_len_p[i]] = '\0';
	  printf("ref: %s (%i)\n", str, input_p->ref_len_p[i]);
	  
	  printf("\n");
     }
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

sw_simd_output_t* sw_simd_output_new(unsigned int depth) {
     sw_simd_output_t* output_p = calloc(1, sizeof(sw_simd_output_t));
     
     output_p->depth = depth;
     output_p->mapped_seq_p = (char**) calloc(depth, sizeof(char*));
     output_p->mapped_ref_p = (char**) calloc(depth, sizeof(char*));
     output_p->mapped_len_p = (unsigned int*) calloc(depth, sizeof(unsigned int));
     //output_p->mapped_ref_len_p = (unsigned int*) calloc(depth, sizeof(unsigned int));
     
     output_p->start_p = (unsigned int*) calloc(depth, sizeof(unsigned int));
     output_p->start_seq_p = (unsigned int*) calloc(depth, sizeof(unsigned int));
     
     output_p->score_p = (float*) calloc(depth, sizeof(float));
     output_p->norm_score_p = (float*) calloc(depth, sizeof(float));
     
     return output_p;
}

//------------------------------------------------------------------------------------

void sw_simd_output_free(sw_simd_output_t* output_p) {
     if (output_p == NULL) return;
     
     if (output_p->mapped_seq_p != NULL) free(output_p->mapped_seq_p);
     if (output_p->mapped_ref_p != NULL) free(output_p->mapped_ref_p);
     if (output_p->mapped_len_p != NULL) free(output_p->mapped_len_p);
     //if (output_p->mapped_ref_len_p != NULL) free(output_p->mapped_ref_len_p);
     
     if (output_p->start_p != NULL) free(output_p->start_p);
     if (output_p->start_seq_p != NULL) free(output_p->start_seq_p);
     if (output_p->score_p != NULL) free(output_p->score_p);
     if (output_p->norm_score_p != NULL) free(output_p->norm_score_p);
     
     free(output_p);
}
//------------------------------------------------------------------------------------

void sw_simd_output_display(unsigned int depth, sw_simd_output_t* output_p) {
     char str[2048];
     
     printf("sw_simd_output (depth = %i of %i):\n", depth, output_p->depth);
     for(int i = 0; i < depth; i++) {
	  
	  memcpy(str, output_p->mapped_seq_p[i], output_p->mapped_len_p[i]);
	  str[output_p->mapped_len_p[i]] = '\0';
	  printf("%i:\n", i);
	  printf("seq: %s (len = %i)\n", str, output_p->mapped_len_p[i]);
	  
	  memcpy(str, output_p->mapped_ref_p[i], output_p->mapped_len_p[i]);
	  str[output_p->mapped_len_p[i]] = '\0';
	  printf("ref: %s (len = %i, start = %i)\n", str, output_p->mapped_len_p[i], output_p->start_p[i]);
	  
	  printf("dif: ");
	  for(int j = 0; j < output_p->mapped_len_p[i]; j++) {
	       printf(output_p->mapped_seq_p[i][j] == output_p->mapped_ref_p[i][j] ? " " : "x");
	  }
	  
	  printf(" (score = %.2f, norm. score = %.2f)\n", output_p->score_p[i], output_p->norm_score_p[i]);
	  
	  printf("\n");
     }
}

//------------------------------------------------------------------------------------

sw_simd_context_t* sw_simd_context_new(float match, float mismatch, float gap_open, float gap_extend) {
  sw_simd_context_t* context_p = (sw_simd_context_t*) malloc(sizeof(sw_simd_context_t));

  context_p->gap_open = gap_open;
  context_p->gap_extend = gap_extend;
  
  context_p->substitution[0] = mismatch;
  context_p->substitution[1] = match;
  
  context_p->x_size = 0;
  context_p->y_size = 0;
  context_p->max_size = 0;
  
  context_p->E = NULL;
  context_p->F = NULL;
  context_p->H = NULL;
  context_p->C = NULL;
  
  context_p->compass_p = NULL;
  
  context_p->a_map = NULL;
  context_p->b_map = NULL;
  
  context_p->q_aux = NULL;
  context_p->r_aux = NULL;
  context_p->aux_size = 0;
  context_p->H_size = 0;
  context_p->F_size = 0; 
  
  // init matrix to -1000.0f
  for (int i = 0; i<128; i++) {
    for (int j = 0; j<128; j++) {
      context_p->matrix[i][j] = -1000.0f;
    }
  }
  
  context_p->matrix['A']['A'] = match;    context_p->matrix['C']['A'] = mismatch; context_p->matrix['T']['A'] = mismatch; context_p->matrix['G']['A'] = mismatch;
  context_p->matrix['A']['C'] = mismatch; context_p->matrix['C']['C'] = match;    context_p->matrix['T']['C'] = mismatch; context_p->matrix['G']['C'] = mismatch;
  context_p->matrix['A']['T'] = mismatch; context_p->matrix['C']['T'] = mismatch; context_p->matrix['T']['T'] = match;    context_p->matrix['G']['T'] = mismatch;
  context_p->matrix['A']['G'] = mismatch; context_p->matrix['C']['G'] = mismatch; context_p->matrix['T']['G'] = mismatch; context_p->matrix['G']['G'] = match;
  context_p->matrix['A']['N'] = mismatch; context_p->matrix['C']['N'] = mismatch; context_p->matrix['T']['N'] = mismatch; context_p->matrix['G']['N'] = mismatch;


  context_p->matrix['N']['A'] = mismatch;
  context_p->matrix['N']['C'] = mismatch;
  context_p->matrix['N']['G'] = mismatch;
  context_p->matrix['N']['N'] = match;
     //     sw_simd_context_update(200, 800, context_p);


     return context_p;
}

//------------------------------------------------------------------------------------

void sw_simd_context_update(int x_size, int y_size, sw_simd_context_t* context_p) {
  //printf("sw_simd_context_update...\n");
     int max_size = x_size * y_size;
     if (context_p->x_size < x_size || context_p->y_size < y_size) {
	  if (context_p->x_size < x_size) {
	       context_p->x_size = x_size;
	  }
	  if (context_p->y_size < y_size) {
	       context_p->y_size = y_size;
	  }
	  context_p->max_size = x_size * y_size;

	  if (context_p->E != NULL) _mm_free(context_p->E);
	  context_p->E = (float*) _mm_malloc(context_p->x_size * SIMD_DEPTH * sizeof(float), SIMD_ALIGN);
	  
	  if (context_p->F != NULL) _mm_free(context_p->F);
	  context_p->F = (float*) _mm_malloc(SIMD_DEPTH * sizeof(float), SIMD_ALIGN);
	  
	  if (context_p->H != NULL) _mm_free(context_p->H);
	  context_p->H = (float*) _mm_malloc(context_p->max_size * SIMD_DEPTH * sizeof(float), SIMD_ALIGN);
	  memset(context_p->H, 0, context_p->max_size * SIMD_DEPTH * sizeof(float));

	  if (context_p->compass_p != NULL) free(context_p->compass_p);
	  context_p->compass_p = (int*) calloc(context_p->max_size * SIMD_DEPTH, sizeof(int));
	  
	  if (context_p->a_map != NULL) free(context_p->a_map);
	  context_p->a_map = (char*) calloc(context_p->max_size + 100, sizeof(char));

	  if (context_p->b_map != NULL) free(context_p->b_map);
	  context_p->b_map = (char*) calloc(context_p->max_size + 100, sizeof(char));
     }
     
     // we must set to zero E and F vectors, however
     // for the H matrix, we only need to set to zero the first row and column
     // and it's done when allocating memory (the first row and column are not
     // overwritten)
     memset(context_p->E, 0, context_p->x_size * SIMD_DEPTH * sizeof(float));
     memset(context_p->F, 0, SIMD_DEPTH * sizeof(float));
     //printf("sw_simd_context_update done !!\n");
}

//------------------------------------------------------------------------------------

void sw_simd_context_free(sw_simd_context_t* context_p) {
  if (context_p == NULL) return;

  if (context_p->E != NULL) _mm_free(context_p->E);  
  if (context_p->H != NULL)  _mm_free(context_p->H); 
  if (context_p->C != NULL) _mm_free(context_p->C);
  if (context_p->F != NULL) _mm_free(context_p->F); 

  if (context_p->compass_p != NULL) free(context_p->compass_p);
  if (context_p->a_map != NULL) free(context_p->a_map);
  if (context_p->b_map != NULL) free(context_p->b_map);
  if (context_p->q_aux != NULL) free(context_p->q_aux);
  if (context_p->r_aux != NULL) free(context_p->r_aux);
  //free(context_p);
  
  free(context_p);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void smith_waterman_simd(sw_simd_input_t* input, sw_simd_output_t* output, 
			 sw_simd_context_t* context) {

  int len, simd_depth = 4, num_queries = input->depth;
  int max_q_len = 0, max_r_len = 0;

  //float gap_open = context->gap_open;
  //float gap_extend = context->gap_extend;

  //  printf("Process %d reads\n", num_queries);


  /*context->H = NULL;
  context->F = NULL;
  context->C = NULL;
  context->H_size = 0;
  context->F_size = 0;*/
  //  context->C_size = 0;

  //printf("Process SW\n");
  for (size_t i = 0; i < num_queries; i++) {
    //printf("REF(%d)(%d)\n", strlen(input->ref_p[i]), strlen(input->seq_p[i]));
    //printf("REF(%d)/(%d) || (%d)/(%d)\n", input->ref_len_p[i], strlen(input->ref_p[i]), input->seq_len_p[i], strlen(input->seq_p[i]));
    //printf("SEQ: %s\n", input->seq_p[i]);
    if (input->seq_len_p[i] > max_q_len){
      max_q_len = input->seq_len_p[i];
      //max_r_len = input->seq_len_p[i];
    }
    if (input->ref_len_p[i] > max_r_len) {
      //max_q_len = input->ref_len_p[i];
      max_r_len = input->ref_len_p[i];
    }
    // if (input->ref_len_p[i] > max_r_len) max_r_len = input->ref_len_p[i];
    //printf("(%i):%s\n", input->seq_len_p[i], input->seq_p[i]);
    //printf("(%i):%s\n", input->ref_len_p[i], input->ref_p[i]);
  }
  

  reallocate_memory(max_q_len, max_r_len, simd_depth, 
		    &context->H_size, &context->H, &context->C, 
		    &context->F_size, &context->F, 
		    &context->aux_size, &context->q_aux, &context->r_aux);

  sse_matrix(num_queries, 
	      input->seq_p, input->seq_len_p, max_q_len, 
	      input->ref_p, input->ref_len_p, max_r_len, 
	      context->matrix, context->gap_open, context->gap_extend, 
	      context->H, context->F, context->C, output->score_p);
  
  simd_traceback(simd_depth, num_queries, 
		  input->seq_p, input->seq_len_p, max_q_len, 
		  input->ref_p, input->ref_len_p, max_r_len, 
		  context->gap_open, context->gap_extend, 
		  context->H, context->C, output->score_p,
		  output->mapped_seq_p, output->start_seq_p,
		  output->mapped_ref_p, output->start_p, 
		  output->mapped_len_p, 
		  context->q_aux, context->r_aux);
  
  // free memory
  /*  printf("Free\n");
  if (context->H != NULL) _mm_free(context->H);
  if (context->C != NULL) _mm_free(context->C);
  if (context->F != NULL) _mm_free(context->F);
  printf("Free end\n");
  */
  for(int i = 0; i < simd_depth; i++){
    //printf("Out: %s\n", output->mapped_seq_p[i]);
    //output->norm_score_p[i] = output->score_p[i] / (input->seq_len_p[i] * context->substitution[1]);
    output->norm_score_p[i] = NORM_SCORE(output->score_p[i], input->seq_len_p[i], context->substitution[1]);
    //printf("Output Score %d, %d, %d\n", output->norm_score_p[i], output->mapped_len_p[i], output->score_p[i]);
  }
}
