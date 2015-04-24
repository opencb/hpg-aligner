#include "sw_server.h"

//====================================================================================
//  Input structure for Smith-Waterman server
//====================================================================================
void sw_optarg_init(float gap_open, float gap_extend, 
		    float match, float mismatch, sw_optarg_t *sw_optarg) {
  sw_optarg->gap_open = gap_open;
  sw_optarg->gap_extend = gap_extend;

  sw_optarg->subst_matrix['A']['A'] = match;
  sw_optarg->subst_matrix['C']['A'] = mismatch;
  sw_optarg->subst_matrix['T']['A'] = mismatch;
  sw_optarg->subst_matrix['G']['A'] = mismatch;
  sw_optarg->subst_matrix['N']['A'] = mismatch;

  sw_optarg->subst_matrix['A']['C'] = mismatch;
  sw_optarg->subst_matrix['C']['C'] = match;
  sw_optarg->subst_matrix['T']['C'] = mismatch;
  sw_optarg->subst_matrix['G']['C'] = mismatch;
  sw_optarg->subst_matrix['N']['C'] = mismatch;

  sw_optarg->subst_matrix['A']['T'] = mismatch;
  sw_optarg->subst_matrix['C']['T'] = mismatch;
  sw_optarg->subst_matrix['T']['T'] = match;
  sw_optarg->subst_matrix['G']['T'] = mismatch;
  sw_optarg->subst_matrix['N']['T'] = mismatch;

  sw_optarg->subst_matrix['A']['G'] = mismatch;
  sw_optarg->subst_matrix['C']['G'] = mismatch;
  sw_optarg->subst_matrix['T']['G'] = mismatch;
  sw_optarg->subst_matrix['G']['G'] = match;
  sw_optarg->subst_matrix['N']['G'] = mismatch;

  sw_optarg->subst_matrix['A']['N'] = mismatch;
  sw_optarg->subst_matrix['C']['N'] = mismatch;
  sw_optarg->subst_matrix['T']['N'] = mismatch;
  sw_optarg->subst_matrix['G']['N'] = mismatch;
  sw_optarg->subst_matrix['N']['N'] = match;

}


void sw_server_input_init(list_t* sw_list, list_t* alignment_list, unsigned int write_size, 
			  float match, float mismatch, float gap_open, float gap_extend, 
			  int min_score, unsigned int flank_length, genome_t* genome, 
			  size_t max_intron_size, int min_intron_size, 
			  size_t seed_max_distance, bwt_optarg_t* bwt_optarg_p, 
			  avls_list_t *avls_list,
			  cal_optarg_t *cal_optarg_p, bwt_index_t *bwt_index_p,
			  metaexons_t *metaexons, linked_list_t *buffer, 
			  linked_list_t *buffer_hc, FILE *f_sa, FILE *f_hc,
			  int pair_mode, sw_server_input_t* input) {
  
  input->sw_list_p = sw_list;
  input->alignment_list_p = alignment_list;
  input->write_size = write_size;
  input->genome_p = genome;
  input->max_intron_size = max_intron_size;
  input->min_intron_size = min_intron_size;
  input->seed_max_distance = seed_max_distance;
  input->bwt_optarg_p =  bwt_optarg_p; 

  // Smith-Waterman parameters
  sw_optarg_init(gap_open, gap_extend, 
		 match, mismatch, &input->sw_optarg);
  input->match = match;
  input->mismatch = mismatch;
  input->gap_open = gap_open;
  input->gap_extend = gap_extend;
  input->min_score = min_score;

  /*
  // this code is done by the previous sw_optarg_init
  input->sw_optarg.gap_open = gap_open;
  input->sw_optarg.gap_extend = gap_extend;

  input->sw_optarg.subst_matrix['A']['A'] = input->match;
  input->sw_optarg.subst_matrix['C']['A'] = input->mismatch;
  input->sw_optarg.subst_matrix['T']['A'] = input->mismatch;
  input->sw_optarg.subst_matrix['G']['A'] = input->mismatch;
  input->sw_optarg.subst_matrix['N']['A'] = input->mismatch;

  input->sw_optarg.subst_matrix['A']['C'] = input->mismatch;
  input->sw_optarg.subst_matrix['C']['C'] = input->match;
  input->sw_optarg.subst_matrix['T']['C'] = input->mismatch;
  input->sw_optarg.subst_matrix['G']['C'] = input->mismatch;
  input->sw_optarg.subst_matrix['N']['C'] = input->mismatch;

  input->sw_optarg.subst_matrix['A']['T'] = input->mismatch;
  input->sw_optarg.subst_matrix['C']['T'] = input->mismatch;
  input->sw_optarg.subst_matrix['T']['T'] = input->match;
  input->sw_optarg.subst_matrix['G']['T'] = input->mismatch;
  input->sw_optarg.subst_matrix['N']['T'] = input->mismatch;


  input->sw_optarg.subst_matrix['A']['G'] = input->mismatch;
  input->sw_optarg.subst_matrix['C']['G'] = input->mismatch;
  input->sw_optarg.subst_matrix['T']['G'] = input->mismatch;
  input->sw_optarg.subst_matrix['G']['G'] = input->match;
  input->sw_optarg.subst_matrix['N']['G'] = input->mismatch;

  input->sw_optarg.subst_matrix['A']['N'] = input->mismatch;
  input->sw_optarg.subst_matrix['C']['N'] = input->mismatch;
  input->sw_optarg.subst_matrix['T']['N'] = input->mismatch;
  input->sw_optarg.subst_matrix['G']['N'] = input->mismatch;
  input->sw_optarg.subst_matrix['N']['N'] = input->match;
  */

  // CAL
  input->flank_length = flank_length;
  input->avls_list = avls_list;

  input->cal_optarg_p = cal_optarg_p;
  input->bwt_index_p = bwt_index_p;
  input->metaexons = metaexons;

  input->buffer = buffer;
  input->buffer_hc = buffer_hc;

  input->f_sa = f_sa;
  input->f_hc = f_hc;

  input->pair_mode = pair_mode;
}

//====================================================================================
//  Smith-Waterman channel for SIMD implementation
//====================================================================================

inline void sw_channel_allocate_ref(unsigned int length, sw_channel_t* channel_p) {
     if (channel_p == NULL) return;
     
     if (channel_p->allocated_ref_size < length) {
	  if (channel_p->ref_p == NULL) {
	       channel_p->ref_p = (char*) calloc(length, sizeof(char));
	  } else {
	       free(channel_p->ref_p);
	       channel_p->ref_p = (char*) calloc(length, sizeof(char));
	  }
	  channel_p->allocated_ref_size = length;
     }
}

//------------------------------------------------------------------------------------

inline void sw_channel_update(size_t read_index, unsigned int cal_index, unsigned int read_len,
			      unsigned int header_len, unsigned int ref_len, sw_channel_t *channel_p) {
     channel_p->read_index = read_index;
     channel_p->cal_index = cal_index;
     channel_p->read_len = read_len;
     channel_p->header_len = header_len;
     channel_p->ref_len = ref_len;
}

//====================================================================================
// main sw function
//====================================================================================

void set_sw_sequences(char **q, char **r, size_t sw_count, char *sequence, 
		      genome_t *genome, int chromosome, seed_region_t *sr) {
  int gap_len;

  // get query sequence, revcomp if necessary
  gap_len = sr->read_end - sr->read_start + 1;
  q[sw_count] = (char *) malloc((gap_len + 1) * sizeof(char));
  memcpy(q[sw_count], &sequence[sr->read_start], gap_len);
  q[sw_count][gap_len] = '\0';

  // get ref. sequence
  gap_len = sr->genome_end - sr->genome_start + 1;
  r[sw_count] = (char *) malloc((gap_len + 1) * sizeof(char));
  genome_read_sequence_by_chr_index(r[sw_count], 0, chromosome, 
  				    &sr->genome_start, &sr->genome_end, genome);
  r[sw_count][gap_len] = '\0';
}

//------------------------------------------------------------------------------------
// apply_sw
//------------------------------------------------------------------------------------

int apply_sw(sw_server_input_t* input, batch_t *batch) {

  mapping_batch_t *mapping_batch = batch->mapping_batch;
  genome_t *genome = input->genome_p;
  sw_optarg_t *sw_optarg = &input->sw_optarg;

  // fill gaps between seeds
  fill_gaps(mapping_batch, sw_optarg, genome, 20, 5);
  merge_seed_regions(mapping_batch);
  fill_end_gaps(mapping_batch, sw_optarg, genome, 20, 400);

  // now we can create the alignments
  fastq_read_t *read;
  array_list_t *fq_batch = mapping_batch->fq_batch;
  
  char *match_seq, *match_qual;
  size_t read_index, read_len, match_len, match_start;
  
  cal_t *cal;
  array_list_t *cal_list = NULL;
  size_t num_cals, num_targets = mapping_batch->num_targets;
  
  seed_region_t *s;
  cigar_code_t *cigar_code;
  cigar_op_t *first_op;

  float score, norm_score, min_score = input->min_score;

  alignment_t *alignment;
  array_list_t *alignment_list;

  char *p, *optional_fields;
  int optional_fields_length, AS;
  
  for (size_t i = 0; i < num_targets; i++) {
    read_index = mapping_batch->targets[i];
    read = (fastq_read_t *) array_list_get(read_index, fq_batch);

    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
    
    if (num_cals <= 0) continue;
    
    read_len = read->length;
    
    alignment_list = array_list_new(num_cals, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {

      // get cal and read index
      cal = array_list_get(j, cal_list);
      if (cal->sr_list->size == 0) continue;

      s = (seed_region_t *) linked_list_get_first(cal->sr_list);
      cigar_code = (cigar_code_t *) s->info;

      norm_score = cigar_code_get_score(read_len, cigar_code);
      score = norm_score * 100; //read_len;
      LOG_DEBUG_F("score = %0.2f\n", norm_score);

      // filter by SW score
      if (norm_score > min_score) {

	// update cigar and sequence and quality strings
	cigar_code_update(cigar_code);
	LOG_DEBUG_F("\tcigar code = %s\n", new_cigar_code_string(cigar_code));
	match_start = 0;
	match_len = cigar_code_nt_length(cigar_code); 
	first_op = cigar_code_get_first_op(cigar_code);
	match_start = (first_op && first_op->name == 'H' ? first_op->number : 0);

	match_seq = (char *) malloc((match_len + 1)* sizeof(char));
	memcpy(match_seq, &read->sequence[match_start], match_len);
	match_seq[match_len] = 0;

	match_qual = (char *) malloc((match_len + 1)* sizeof(char));
	memcpy(match_qual, &read->quality[match_start], match_len);
	match_qual[match_len] = 0;

	// set optional fields
	optional_fields_length = 100;
	optional_fields = (char *) calloc(optional_fields_length, sizeof(char));
      
	p = optional_fields;
	AS = (int) norm_score * 100;
	
	sprintf(p, "ASi");
	p += 3;
	memcpy(p, &AS, sizeof(int));
	p += sizeof(int);
	
	sprintf(p, "NHi");
	p += 3;
	memcpy(p, &num_cals, sizeof(int));
	p += sizeof(int);
	
	sprintf(p, "NMi");
	p += 3;
	memcpy(p, &cigar_code->distance, sizeof(int));
	p += sizeof(int);

	assert(strlen(read->sequence) == cigar_code_nt_length(cigar_code));

	// create an alignment and insert it into the list
	alignment = alignment_new();

	//read_id = malloc(read->length);
	size_t header_len = strlen(read->id);
	char *head_id = (char *) malloc(header_len + 1);

	get_to_first_blank(read->id, header_len, head_id);

	alignment_init_single_end(head_id, match_seq, match_qual, 
				  cal->strand, cal->chromosome_id - 1, cal->start - 1,
				  new_cigar_code_string(cigar_code), 
				  cigar_code_get_num_ops(cigar_code), 
				  norm_score * 254, 1, (num_cals > 1),
				  optional_fields_length, (uint8_t *)optional_fields, alignment);

	array_list_insert(alignment, alignment_list);
      }
    }
    
    // free the cal list, and update the mapping list with the alignment list
    array_list_free(cal_list, (void *) cal_free);
    mapping_batch->mapping_lists[read_index] = alignment_list;
  }
  // go to the next stage
  return DNA_POST_PAIR_STAGE;
}

//--------------------------------------------------------------------------------------


void fill_matrix(subst_matrix_t subst_matrix, float match, float mismatch, int type, float factor_match, float factor_mismatch) {

  subst_matrix['A']['A'] = match;
  subst_matrix['C']['A'] = mismatch;
  subst_matrix['T']['A'] = mismatch;
  subst_matrix['N']['A'] = mismatch;

  subst_matrix['A']['C'] = mismatch;
  subst_matrix['G']['C'] = mismatch;
  subst_matrix['N']['C'] = mismatch;

  subst_matrix['A']['T'] = mismatch;
  subst_matrix['T']['T'] = match;
  subst_matrix['G']['T'] = mismatch;
  subst_matrix['N']['T'] = mismatch;

  subst_matrix['C']['G'] = mismatch;
  subst_matrix['T']['G'] = mismatch;
  subst_matrix['N']['G'] = mismatch;

  subst_matrix['A']['N'] = mismatch;
  subst_matrix['C']['N'] = mismatch;
  subst_matrix['T']['N'] = mismatch;
  subst_matrix['G']['N'] = mismatch;
  subst_matrix['N']['N'] = match;

  /*
  float a = match * factor_match;
  float b = mismatch / factor_mismatch;
  float c = mismatch / factor_mismatch;
  */
  float x = 5;
  float y = 5;
  float z = 5;

  if (type == 1) {
    subst_matrix['C']['C'] = x;
    subst_matrix['T']['C'] = y;
    subst_matrix['C']['T'] = z;
    subst_matrix['G']['G'] = match;
    subst_matrix['A']['G'] = mismatch;
    subst_matrix['G']['A'] = mismatch;
  } else {
    subst_matrix['C']['C'] = match;
    subst_matrix['T']['C'] = mismatch;
    subst_matrix['C']['T'] = mismatch;
    subst_matrix['G']['G'] = x;
    subst_matrix['A']['G'] = y;
    subst_matrix['G']['A'] = z;
  }

}


//--------------------------------------------------------------------------------------

sw_output_t *sw_output_new(int strand, size_t chrom, size_t ref_start, size_t ref_len,
			   size_t mref_len, size_t mquery_start, size_t mref_start,
			   float score, float norm_score, char* mquery, char* mref) {

  sw_output_t *p = (sw_output_t *) calloc(1, sizeof(sw_output_t));

  p->strand = strand;
  p->chromosome = chrom;
  p->ref_start = ref_start;
  p->ref_len = ref_len;
  p->mref_len = mref_len;
  p->mquery_start = mquery_start;
  p->mref_start = mref_start;
  p->score = score;
  p->norm_score = norm_score;
  p->mquery = strdup(mquery);
  p->mref = strdup(mref);
  //p->mquery = NULL;
  //p->mref = NULL;

  return p;
}

//--------------------------------------------------------------------------------------

void sw_output_free(sw_output_t *p) {
  if (p == NULL) return;

  if (p->mquery != NULL) free(p->mquery);
  if (p->mref != NULL) free(p->mref);

  free(p);
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

