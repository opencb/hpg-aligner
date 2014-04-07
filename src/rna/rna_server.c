#include <omp.h>
#include "rna_server.h"

//COLORS
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

#define CALING_DOUBLE_ANCHORS  0
#define CALING_LEFT_ANCHORS    1
#define CALING_RIGHT_ANCHORS   2

#define META_OPEN  0
#define META_CLOSE 1

#define CIGAR_SW_MIDDLE      0
#define CIGAR_SIMPLE_MIDDLE  1
#define CIGAR_ANCHOR_LEFT    2
#define CIGAR_ANCHOR_RIGHT   3

#define MAX_DEPTH 4

#define SIMPLE_SW       1
#define SP_SW           2
#define EXTREM_SW_LEFT  4
#define EXTREM_SW_RIGHT 5

#define SP_METAEXON     3

#define SW_NORMAL 0
#define SW_FINAL 1

#define OP_TYPE 0
#define REFERENCE_LEFT 1
#define REFERENCE_MIDLE 2
#define REFERENCE_INTRON 3
#define REFERENCE_RIGHT 4

#define MINIMUN_CAL_LENGTH 10
#define ERRORS_ZONE 8
#define MINIMUM_INTRON_LENGTH 10
#define MAX_CIGAR_OPERATIONS 20

#define LEFT_EXTREM 0
#define RIGHT_EXTREM 1

//===============================================

#define A_NT  0
#define C_NT  1
#define G_NT  2
#define T_NT  3

const unsigned char splice_nt[4] = {'A', 'C', 'G', 'T'};

//--------------------------------------------//
//            Not found splice junction       //
//--------------------------------------------//

#define NOT_SPLICE	-1

//--------------------------------------------//
//      No Cannonical Splice junction         //
//--------------------------------------------//

#define UNKNOWN_SPLICE	0
                                            
//--------------------------------------------//
//        Cannonical Splice Junction          //
//--------------------------------------------//

#define GT_AG_SPLICE  	1 //+
#define CT_AC_SPLICE  	2 //-
  
//--------------------------------------------//
//      Semi-Cannonical Splice Junction       //
//--------------------------------------------//

#define AT_AC_SPLICE  	3 //+
#define GT_AT_SPLICE  	4 //-
#define GC_AG_SPLICE  	5 //+
#define CT_GC_SPLICE  	6 //-

//===============================================


#ifndef MAX
   #define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

//======== SPLICE JUNCTION TYPE ============
//#define NOT_SPLICE	        0
//#define GT_AG_SPLICE  	1
//#define CT_AC_SPLICE  	2
//==========================================

#define SEARCH_START 0
#define SEARCH_END 1
#define POSSIBLES_MARKS 2
#define STRANDS_NUMBER 2
#define MAX_FUSION 2028
#define MIN_HARD_CLIPPING 10


#define SW_LEFT   0
#define SW_RIGHT  1
#define SW_MIDDLE 2

int ii = -1;

//extern size_t TOTAL_READS_PROCESS, TOTAL_SW, TOTAL_READS_SA;
extern pthread_mutex_t mutex_sp; 

//extern size_t tot_reads_in;
//extern size_t tot_reads_out;


cigar_code_t *search_left_single_anchor(int gap_close, 
					cal_t *cal,
					int filter_pos, 
					array_list_t *right_breaks,
					char *query_map,
					metaexons_t *metaexons,
					genome_t *genome, 
					avls_list_t *avls_list);

cigar_code_t *search_right_single_anchor(int gap_close, 
					 cal_t *cal,
					 int filter_pos, 
					 array_list_t *left_breaks,
					 char *query_map, metaexons_t *metaexons,
					 genome_t *genome, 
					 avls_list_t *avls_list);

char cigar_automata_status(unsigned char status) {

  switch (status) {
    case CIGAR_MATCH_MISMATCH:
      return 'M';
      break;
    case CIGAR_INSERTION:
      return 'I';
      break;
    case CIGAR_DELETION:
      return 'D';
      break;
    case CIGAR_SKIPPED:
      return 'N';
      break;
    case CIGAR_PADDING:
      return 'P';
      break;
  }

  return ' ';
}

// Select the splice junction type
int splice_junction_type(char nt_start_1, char nt_start_2, char nt_end_1, char nt_end_2) {
  //LOG_DEBUG_F("SEARCH SPLICE JUNCTION TYPE FOR %c%c - %c%c\n", nt_start_1, nt_start_2, nt_end_1, nt_end_2);

  int splice_type = NOT_SPLICE;

  if (nt_start_1 == splice_nt[G_NT]) {
    if (nt_start_2 == splice_nt[T_NT]) {
      if (nt_end_1 == splice_nt[A_NT]) {
	if (nt_end_2 == splice_nt[G_NT]) {
	  //Report GT-AG Splice
	  splice_type = GT_AG_SPLICE;
	} else if (nt_end_2 == splice_nt[T_NT]) {
	  //Report GT-AT Splice
	  splice_type = GT_AT_SPLICE;
	}
      }
    } else if (nt_start_2 == splice_nt[C_NT] && 
	       nt_end_1 == splice_nt[A_NT] && 
	       nt_end_2 == splice_nt[G_NT]) {
      //Report GC -AG
      splice_type = GC_AG_SPLICE;
    }
  } else if (nt_start_1 == splice_nt[C_NT] && 
	     nt_start_2 == splice_nt[T_NT]) {
    if (nt_end_1 == splice_nt[A_NT] && 
	nt_end_2 == splice_nt[C_NT]) {
      //Report CT-AC
      splice_type = CT_AC_SPLICE;
    }else if (nt_end_1 == splice_nt[G_NT] && 
	      nt_end_2 == splice_nt[C_NT]) {
      //Report CT-GC
      splice_type = CT_GC_SPLICE;
    } 
  } else if (nt_start_1 == splice_nt[A_NT] && 
	     nt_start_2 == splice_nt[T_NT] && 
	     nt_end_1 == splice_nt[A_NT] && 
	     nt_end_2 == splice_nt[C_NT] ) {
    //Report AT-AC
    splice_type = AT_AC_SPLICE;
  }

  return splice_type;

}


cigar_code_t *generate_cigar_sw_output(char *seq_sw, 
				       char *ref_sw,
				       size_t l_exon_start,
				       size_t l_exon_end,
				       size_t r_exon_start,
				       size_t r_exon_end,
				       int chromosome,
				       int strand,
				       int seq_start, 
				       int ref_start,
				       int len_orig_seq,
				       int len_orig_ref,
				       avls_list_t *avls_list,
				       avl_node_t **node_avl_start,
				       avl_node_t **node_avl_end, 
				       genome_t *genome,
				       int min_intron_size) {

  //printf("[%lu-%lu] - [%lu-%lu]\n", l_exon_start, l_exon_end, r_exon_start, r_exon_end);
  //printf("Ref: %s\n", ref_sw);
  //printf("Seq: %s\n", seq_sw);

  *node_avl_start = NULL;
  *node_avl_end   = NULL;

  int const MIN_GAP_SEARCH = 15;
  int const EXTRA_SEARCH = 5;
  int map_sw_len = strlen(seq_sw);
  unsigned char automata_status = CIGAR_MATCH_MISMATCH;
  int found = NOT_SPLICE;
  int cigar_value = 0;
  int j = 0;
  int start_gap, end_gap;
  cigar_op_t *cigar_op;  
  int insertions_tot = 0;
  int deletions_tot = 0;
  int tot_matches = 0;
  int padding_left = ref_start;
  //int mode = MODE_EXON;
  int len_ref, len_r_gap;
  int gap_len;
  int cnt_ext = 0;
  int pos;
  int sw_seq_len = 0;
  
  size_t start_splice, end_splice;
  cigar_code_t *cigar_code = cigar_code_new();
  cigar_op_t *op;
  
  int n_splice = 0;
  
  if (seq_start > 0) {
    //Middle or last ref
    if (ref_start == 0) {
      cigar_code_append_op(cigar_op_new(seq_start, 'I'), cigar_code);
    } else {
      if (ref_start == seq_start) {
	cigar_code_append_op(cigar_op_new(seq_start, 'M'), cigar_code);
      } else {
	if (ref_start > seq_start) {
	  cigar_code_append_op(cigar_op_new(ref_start - seq_start, 'D'), cigar_code);
	  cigar_code_append_op(cigar_op_new(seq_start, 'M'), cigar_code);
	  cigar_code->distance += seq_start;
	} else {
	  cigar_code_append_op(cigar_op_new(seq_start - ref_start, 'I'), cigar_code);
	  cigar_code_append_op(cigar_op_new(ref_start, 'M'), cigar_code);
	  cigar_code->distance += ref_start;
	} 
      }
    }
  } else if (ref_start > 0) {
    cigar_code_append_op(cigar_op_new(ref_start, 'D'), cigar_code);
  }

  while (j < map_sw_len) {
    //printf("NT-SEQ (%c.%c)\n", seq_sw[j], ref_sw[j]);
    if (ref_sw[j] != '-'  && seq_sw[j] != '-') {
      tot_matches++;
      padding_left++;
      //printf("\tpadding left = %i (%c-%c)\n", padding_left, ref_sw[j], seq_sw[j]);
      if (ref_sw[j] != seq_sw[j]) {
	cigar_code->distance++;
      }
      //Match/Mismatch Area	  
      if (automata_status == CIGAR_MATCH_MISMATCH) {
	cigar_value++;
	//printf("\t++%i (%c-%c)\n", cigar_value, ref_sw[j], seq_sw[j]);
      } else {
	//printf("Report %i\n", cigar_value);
	cigar_code_append_new_op(cigar_value, cigar_automata_status(automata_status), cigar_code);
	automata_status = CIGAR_MATCH_MISMATCH;
	cigar_value = 1;
      } 
    } else if (ref_sw[j] == '-' && seq_sw[j] != '-') {
      insertions_tot++;
      //Insertion Area
      if (automata_status == CIGAR_INSERTION) {
	cigar_value++;
      } else {
	//printf("Report %i\n", cigar_value);
	cigar_code_append_new_op(cigar_value, cigar_automata_status(automata_status), cigar_code);
	automata_status = CIGAR_INSERTION;
	cigar_value = 1;
      }
    } else if (ref_sw[j] != '-' && seq_sw[j] == '-') {
      //printf("Report %i\n", cigar_value);
      //op = cigar_op_new(cigar_value, cigar_automata_status(automata_status));
      cigar_op_t *op_aux = array_list_get(cigar_code->ops->size - 1, cigar_code->ops);
      if (op_aux != NULL && op_aux->name == cigar_automata_status(automata_status)) { 
	op_aux->number += cigar_value;//op->number;
	op = op_aux;
      } else {
	op = cigar_op_new(cigar_value, cigar_automata_status(automata_status));
	cigar_code_append_op(op, cigar_code);
      }
      
      //Deletion Area. Travel in the deletions gap to found some splice junction
      start_gap = j;
      while (j < map_sw_len && seq_sw[j] == '-' ) {	  
	j++;
	deletions_tot++;
      }
            
      end_gap = j - 1;
      gap_len = end_gap - start_gap + 1;
      //printf("gap_len = %i\n", gap_len);

      //Search gap start and gap end
      if (gap_len > MIN_GAP_SEARCH) {
	if (j >= map_sw_len || start_gap + 1 >= map_sw_len) {
	  //printf("ERROR OVERFLOW SW id.1");
	  array_list_clear(cigar_code->ops, (void *)cigar_op_free);
	  cigar_code_free(cigar_code); 
	  return NULL; 
	}
	//Search splice junction
	found = splice_junction_type(ref_sw[start_gap], ref_sw[start_gap + 1], ref_sw[end_gap - 1], ref_sw[end_gap]);
	//printf("Found %i == %i\n", found, NOT_SPLICE);

	if (found == NOT_SPLICE) {
	  //Search Xnt (+)---->	
	  cnt_ext = 1;
	  while (cnt_ext < EXTRA_SEARCH) {
	    if (start_gap + cnt_ext + 1 >= map_sw_len || 
		end_gap + cnt_ext >= map_sw_len) {
	      //printf("ERROR OVERFLOW SW id.2");
	      array_list_clear(cigar_code->ops, (void *)cigar_op_free);
	      cigar_code_free(cigar_code); 
	      return NULL;
	    }
	    found = splice_junction_type(ref_sw[start_gap + cnt_ext], ref_sw[start_gap + cnt_ext + 1], 
					 ref_sw[end_gap + cnt_ext - 1], ref_sw[end_gap + cnt_ext]);	       
	    if (found != NOT_SPLICE) {
	      //printf("Found %i\n", found);
	      break;
	    } else {
	      //printf("continue search...\n");
	    }
	    cnt_ext++;
	  }
	}
	
	//if (!found) { printf("Seq: %s\n", seq_sw); printf("Ref: %s\n", ref_sw); }
	char str_sp_type[10];
	if (found == NOT_SPLICE) {
	  cnt_ext = 0;
	}

	len_ref = l_exon_end - l_exon_start + 1;
	len_r_gap = gap_len - (len_ref - padding_left/*(tot_matches + ref_start)*/);
	if (len_r_gap < 0) { assert(len_r_gap); }
	LOG_DEBUG_F("Calculating SP: l_exon_end = %lu, l_exon_start = %lu, gap_len = %i, len_ref = %i, tot_matches = %i, len_r_gap=%i, padding_left = %i, cnt_ext = %i, seq_start = %i, ref_start = %i, r_exon_start = %i\n",
		    l_exon_end, l_exon_start, gap_len, len_ref, tot_matches, len_r_gap, padding_left, cnt_ext, seq_start, ref_start, r_exon_start);
	
	//printf("cnt_ext=%i\n", cnt_ext);
	start_splice = l_exon_start + padding_left + cnt_ext;
	end_splice = r_exon_start + len_r_gap - 1 + cnt_ext;
	cigar_value = end_splice - start_splice;

	//printf("%i - %i = %i\n", start_splice, end_splice, cigar_value);

	if (cigar_value < 0) { return NULL; }
	
	cigar_code_append_new_op(cigar_value + 1, 'N', cigar_code);	
	if (found == CT_AC_SPLICE || found == GT_AT_SPLICE || found == CT_GC_SPLICE ) {
	  strand = 1;
	} else if (found != NOT_SPLICE) { 
	  strand = 0;
	} else {
	  /*char nt_start[5], nt_end[5];
	  size_t g_start = start_splice;
	  size_t g_end   = start_splice + 1;
	  genome_read_sequence_by_chr_index(nt_start, strand, chromosome - 1, &g_start, &g_end, genome);
	  
	  g_end   = end_splice;
	  g_start = end_splice - 1;
	  genome_read_sequence_by_chr_index(nt_end, strand, chromosome - 1, &g_start, &g_end, genome);
	  
	  found = UNKNOWN_SPLICE;
	  sprintf(str_sp_type, "%c%c-%c%c\0", nt_start[0], nt_start[1], nt_end[0], nt_end[1]);
	  */
	  //printf("SW:(%i)[%i:%lu-%lu] : %c%c vs %c%c\n", strand, chromosome, start_splice, end_splice, nt_start[0], nt_start[1], nt_end[0], nt_end[1]);
	}

	//printf("SP COORDS(%i)=> [%i:%lu-%lu] = %d\n", strand, chromosome, start_splice, end_splice, cigar_value);

	//if (chromosome == 1 && start_splice == 17743 && end_splice == 17914) { printf("%s ::-->\n", id); exit(-1); }
	n_splice++;
	if (n_splice > 1) { 
	  array_list_clear(cigar_code->ops, (void *)cigar_op_free);
	  cigar_code_free(cigar_code); 
	  return NULL; 
	}

	if (found != NOT_SPLICE) {
	  if (end_splice - start_splice + 1 >= min_intron_size) {
	    //printf("SP:=>%i:%lu-%lu\n", chromosome, start_splice, end_splice);
	    allocate_start_node(chromosome - 1,
				strand,
				start_splice,
				end_splice,
				start_splice,
				end_splice,
				FROM_READ,
				found,
				str_sp_type, 
				node_avl_start,
				node_avl_end, 
				avls_list);
	  } else {
	    array_list_clear(cigar_code->ops, (void *)cigar_op_free);
	    cigar_code_free(cigar_code); 
	    return NULL; 
	  }
	}

	//printf("End allocate %i -> %i\n", op->number, op->number + cnt_ext);
	op->number += cnt_ext;
	  //} else {
	  //printf(":( NOT FOUND! [%lu-%lu]---[%lu-%lu]\n", l_exon_start, l_exon_end, r_exon_start, r_exon_end);
	  //cigar_value = gap_len;	
	  //cigar_code_append_new_op(cigar_value, 'D', cigar_code);	
	  //padding_left += cigar_value;
	  //}
      } else { //Deletions Section
	cigar_value = gap_len;	
	//printf("Report %i\n", cigar_value);
	cigar_code_append_new_op(cigar_value, 'D', cigar_code);	
	
	//printf("padding left = %i + %i\n", padding_left, cigar_value);
	padding_left += cigar_value;
	
	//j--; // ==NEW!== go to prev position
      }
      
      if (j >= map_sw_len) { 
	array_list_clear(cigar_code->ops, (void *)cigar_op_free);
	cigar_code_free(cigar_code);
	return NULL; 
      }

      if (ref_sw[j] != '-') { automata_status = CIGAR_MATCH_MISMATCH; tot_matches++; }
      else { automata_status = CIGAR_INSERTION; insertions_tot++; }
      cigar_value = 1 - cnt_ext;
      cnt_ext = 0;

    } else {
      //Padding Area
      insertions_tot++;
      deletions_tot++;
      if (automata_status == CIGAR_PADDING) {
	cigar_value++;
      } else {
	//printf("Report %i\n", cigar_value);
	cigar_code_append_new_op(cigar_value, cigar_automata_status(automata_status), cigar_code);
	automata_status = CIGAR_PADDING;
	cigar_value = 1;
      }
    }
    j++;
  }

  cigar_code_append_new_op(cigar_value, cigar_automata_status(automata_status), cigar_code);

  //sw_seq_len = tot_matches + tot_insertions + seq_start;

  int map_seq_len  = ((map_sw_len - deletions_tot) + seq_start);
  int map_ref_len  = ((map_sw_len - insertions_tot) + ref_start);
  int last_h, last_h_aux;

  //printf("query_start = %i, ref_start = %i, map_seq_len = %i, map_ref_len = %i, query_len = %i, ref_len = %i, map_len = %i\n", 
  //	 query_start, ref_start, map_seq_len, map_ref_len, query_len, ref_len, map_len);
  //printf("***Map_seq_len = %i, map_Ref_len = %i, len_orig_seq = %i, len_orig_ref = %i\n", 
  //	 map_seq_len, map_ref_len, len_orig_seq, len_orig_ref);

  if (map_seq_len < len_orig_seq) {
    last_h = len_orig_seq - map_seq_len;
    //printf("last_h = %i\n", last_h);
    //Middle or first ref
    if (map_ref_len == len_orig_ref) {
      cigar_code_append_op(cigar_op_new(last_h, 'I'), cigar_code);
    } else {
      last_h_aux = len_orig_ref - map_ref_len;
      //printf("last_h_aux = %i\n", last_h_aux);
      if (last_h_aux == last_h) {
	cigar_code_append_op(cigar_op_new(last_h, 'M'), cigar_code);
      } else {	  
	if (last_h_aux > last_h) {
	  cigar_code_append_op(cigar_op_new(last_h_aux - last_h, 'D'), cigar_code);
	  cigar_code_append_op(cigar_op_new(last_h, 'M'), cigar_code);
	  cigar_code->distance += last_h;
	} else {
	  cigar_code_append_op(cigar_op_new(last_h - last_h_aux, 'I'), cigar_code);
	  cigar_code_append_op(cigar_op_new(last_h_aux, 'M'), cigar_code);
	  cigar_code->distance += last_h_aux;
	} 
      }
    }  
  } else if (map_ref_len < len_orig_ref) {
    cigar_code_append_op(cigar_op_new(len_orig_ref - map_ref_len, 'D'), cigar_code);
  }

  if (n_splice == 0) {
    array_list_clear(cigar_code->ops, (void *)cigar_op_free);
    cigar_code_free(cigar_code);
    return NULL;
  }

  //Refresh distance cigar
  for (int i = 0; i < cigar_code_get_num_ops(cigar_code); i++) {
    cigar_op_t *cigar_op = array_list_get(i, cigar_code->ops);
    if (cigar_op->name == 'D' || cigar_op->name == 'I') {
      cigar_code->distance += cigar_op->number;
    }
  }
  
  //printf(" :::: => %s\n", new_cigar_code_string(cigar_code));

  return cigar_code;

}

typedef struct fusion_coords {
  size_t l_exon_start;
  size_t l_exon_end;
  size_t r_exon_start;
  size_t r_exon_end;
  size_t read_start;
  size_t read_end;
  int chromosome;
  int strand;
  int type_sw;
  size_t read_id;
  char *id;
  cal_t *cal_ref;
  int cal_group_id;
} fusion_coords_t;

fusion_coords_t *fusion_coords_new(size_t l_exon_start,
				   size_t l_exon_end,
				   size_t r_exon_start,
				   size_t r_exon_end,
				   size_t read_start,
				   size_t read_end, 
				   int chromosome,
				   int strand,
				   int type_sw,
				   char *id, 
				   cal_t *cal_ref,
				   size_t read_id,
				   int cal_group_id) {

  fusion_coords_t *fusion_coords = (fusion_coords_t *)malloc(sizeof(fusion_coords_t));

  fusion_coords->l_exon_start = l_exon_start;
  fusion_coords->l_exon_end   = l_exon_end;
  fusion_coords->r_exon_start = r_exon_start;
  fusion_coords->r_exon_end   = r_exon_end;
  fusion_coords->read_start   = read_start;
  fusion_coords->read_end     = read_end;
  fusion_coords->chromosome   = chromosome;
  fusion_coords->strand       = strand;
  fusion_coords->type_sw      = type_sw;
  fusion_coords->id           = id;
  fusion_coords->cal_ref      = cal_ref;
  fusion_coords->read_id      = read_id;
  fusion_coords->cal_group_id = cal_group_id;

  return fusion_coords;
}

void fusion_coords_free(fusion_coords_t *fusion_coords) {
  free(fusion_coords);
}

    //array_list_size(cals_targets[target])
    /*for (i = 0; i < number_of_best; i++) {
      fusion_cals = array_list_get(i, cals_targets[target]);
      printf("\t FUSION PACKAGE %i(%f):\n", i, cals_score[i]);
      for (j = 0; j < array_list_size(fusion_cals); j++) {
	cal = array_list_get(j, fusion_cals);
	printf("\t\t CAL(%i) [%i:%lu-%lu]:\n", j, cal->chromosome_id, cal->start,
	       cal->end);
	
	linked_list_iterator_init(cal->sr_list, &itr);
	s = (seed_region_t *) linked_list_iterator_curr(&itr);
	while (s != NULL) {
	  printf("\t\t\t seed: [%i|%i - %i|%i]\n",
		 s->genome_start, s->read_start, s->read_end, s->genome_end);
	  s = (seed_region_t *) linked_list_iterator_next(&itr);
	}
      } 
      }*/

void extend_by_mismatches(char *ref_e1, char *ref_e2, char *query, 
			  int r_e1, int r_e2, 
			  int q_e1, int q_e2, int lim_err,
			  int *dsp_e1, int *dsp_e2) {
  int len_ref_e1 = strlen(ref_e1);
  int len_ref_e2 = strlen(ref_e2);
  int len_query = strlen(query);
  int num_err = 0;
  int tot_err = 0;
  *dsp_e1 = 0;
  *dsp_e2 = 0;
  //printf("e1=%i, %i, %i, %c == %c, %i, %i, %i, %i, %i\n", e1, strlen(query), 
  //strlen(reference_prev), reference_prev[e1], query[e1], genome_start, genome_end, read_start, read_end, seeds_nt);
  //printf("Ref Ex1: %s\n", &ref_e1[0]);
  //printf("Ref Ex2: %s\n", &ref_e2[0]);
  //printf("Query  : %s\n", &query[0]);
 
  while (q_e1 < len_query && r_e1 < len_ref_e1) {
    //printf("START: %c != %c (%i)\n", ref_e1[r_e1], query[q_e1], num_err);
    if (ref_e1[r_e1] != query[q_e1]) {
      num_err++;
      if (num_err >= lim_err) { break; }
    }
    r_e1++;
    q_e1++;
    (*dsp_e1)++;
  }

  tot_err += num_err;
  num_err = 0;

  while (r_e2 >= 0 && q_e2 >= 0) { 
    //printf("END: %c != %c(%i)\n", ref_e2[r_e2], query[q_e2], num_err);
    if (ref_e2[r_e2] != query[q_e2]) {
      num_err++;
      if (num_err >= lim_err) { break; }
    }
    r_e2--;
    q_e2--;
    (*dsp_e2)++;
  }

  tot_err += num_err;

}

cigar_code_t *meta_alignment_fill_extrem_gap(char *query, 
					     cal_t *cal,
					     int type,
					     genome_t *genome,
					     metaexons_t *metaexons, 
					     avls_list_t *avls_list) {

  //return NULL;
  //printf("FILL EXTREM GAPS...%s\n", type == FILL_GAP_LEFT? "LEFT" : "RIGHT");
  int chromosome_id = cal->chromosome_id;
  size_t genome_start, genome_end;
  int read_start, read_end;
  int read_gap = 0;
  char reference[2048];
  int type_search;
  cigar_code_t *cigar_code = NULL;
  int length = strlen(query);
  metaexon_t *metaexon;

  //cal_print(cal);

  pthread_mutex_lock(&metaexons->mutex[cal->chromosome_id - 1]);

  //FILL_GAP_LEFT  0
  //FILL_GAP_RIGHT 1
  if (type == FILL_GAP_LEFT) {
    //printf("FILL_GAP_LEFT\n");
    seed_region_t *s_prev = linked_list_get_first(cal->sr_list);
    assert(s_prev != NULL);
    read_start = 0;
    read_end = s_prev->read_start - 1;
    read_gap = s_prev->read_start;

    if (read_gap >= s_prev->genome_start) {
      genome_start = 0;
    } else {
      genome_start = s_prev->genome_start - read_gap;
    }

    genome_end = s_prev->genome_start - 1;
    //printf("FILL_GAP_LEFT [%i-%i] %i\n", read_start, read_end, read_gap);
    
    metaexon_search(cal->strand, cal->chromosome_id - 1,
		    cal->start, cal->start + 5, &metaexon,
		    metaexons);
    
  } else {
    seed_region_t *s_prev = linked_list_get_last(cal->sr_list);
    assert(s_prev != NULL);
    read_start = s_prev->read_end + 1;
    read_end = length - 1;
    read_gap = length - s_prev->read_end - 1;
    genome_start = s_prev->genome_end + 1;
    genome_end = s_prev->genome_end + read_gap + 1;
    //printf("FILL_GAP_RIGHT [%i-%i] %i\n", read_start, read_end, read_gap);
    metaexon_search(cal->strand, cal->chromosome_id - 1,
		    cal->end - 5, cal->end, &metaexon,
		    metaexons);
    
  }

  //printf("FILL EXTREM GAP: [%i-%i]\n", read_start, read_end);

  //printf("%i:%lu-%lu(%i)\n", cal->chromosome_id, cal->start, cal->end, cal->strand);

  if (metaexon != NULL) {
    //printf("METAEXON NOT NULL! %i-%i\n", metaexon->start, metaexon->end);
    if (type == FILL_GAP_LEFT) {
      type_search = 2;
    } else {
      type_search = 3;
    }
    //if (type_search == 2 &&
    //!metaexon->left_closed) {
    //type_search = 0;
    //} else if (type_search == 3 &&
    //	       !metaexon->right_closed) {
    //type_search = 0;
    //}
  } else {
    type_search = 0;
  }
  
  
  if (type_search == 2) {
    //Search other
    //printf("Search WITH RIGHT ANCHOR\n");
    cigar_code = search_right_single_anchor(read_gap, 
					    cal,
					    0, 
					    metaexon->left_breaks,
					    query,
					    metaexons,
					    genome, 
					    avls_list);    
  } else if (type_search == 3) {
    //printf("Search WITH LEFT ANCHOR\n");
    cigar_code = search_left_single_anchor(read_gap, 
					   cal,
					   0, 
					   metaexon->right_breaks,
					   query,
					   metaexons,
					   genome, 
					   avls_list);   
  }
  
  pthread_mutex_unlock(&metaexons->mutex[cal->chromosome_id - 1]);


  return cigar_code;

}


cigar_code_t *fill_extrem_gap(char *query, 
			      cal_t *cal,
			      int type,
			      genome_t *genome,
			      metaexons_t *metaexons, 
			      avls_list_t *avls_list) {

  //return NULL;
  //printf("FILL EXTREM GAPS...%s\n", type == FILL_GAP_LEFT? "LEFT" : "RIGHT");
  int chromosome_id = cal->chromosome_id;
  size_t genome_start, genome_end;
  int read_start, read_end;
  int read_gap = 0;
  char reference[2048];
  int type_search;
  cigar_code_t *cigar_code = NULL;
  int length = strlen(query);
  metaexon_t *metaexon;

  //FILL_GAP_LEFT  0
  //FILL_GAP_RIGHT 1
  if (type == FILL_GAP_LEFT) {
    //printf("FILL_GAP_LEFT\n");
    seed_region_t *s_prev = linked_list_get_first(cal->sr_list);
    assert(s_prev != NULL);
    read_start = 0;
    read_end = s_prev->read_start - 1;
    read_gap = s_prev->read_start;

    if (read_gap >= s_prev->genome_start) {
      genome_start = 0;
    } else {
      genome_start = s_prev->genome_start - read_gap;
    }

    genome_end = s_prev->genome_start - 1;
    //printf("FILL_GAP_LEFT [%i-%i] %i\n", read_start, read_end, read_gap);
  } else {
    seed_region_t *s_prev = linked_list_get_last(cal->sr_list);
    assert(s_prev != NULL);
    read_start = s_prev->read_end + 1;
    read_end = length - 1;
    read_gap = length - s_prev->read_end - 1;
    genome_start = s_prev->genome_end + 1;
    genome_end = s_prev->genome_end + read_gap + 1;
    //printf("FILL_GAP_RIGHT [%i-%i] %i\n", read_start, read_end, read_gap);
  }

  //printf("FILL EXTREM GAP: [%i-%i]\n", read_start, read_end);

  pthread_mutex_lock(&metaexons->mutex[cal->chromosome_id - 1]);

  //printf("genome_start = %i, genome_end = %i\n", genome_start, genome_end);

  //printf("%i:%lu-%lu(%i)\n", cal->chromosome_id, cal->start, cal->end, cal->strand);

  metaexon_search(cal->strand, cal->chromosome_id - 1,
		  cal->start, cal->end, &metaexon,
		  metaexons);

  if (metaexon != NULL) {
    //printf("METAEXON NOT NULL! %i-%i\n", metaexon->start, metaexon->end);
    if (type == FILL_GAP_LEFT) {
      //printf("genome_start(%i) >= metaexon->start(%i)\n", genome_start, metaexon->start);
      //if (genome_start >= metaexon->start) {
      //type_search = 1;
      //} else { type_search = 2; }
      type_search = 2;
    } else {
      //printf("genome_end(%i) >= metaexon->end(%i)\n", genome_end, metaexon->end);
      //if (genome_end <= metaexon->end) {
      //type_search = 1;
      //} else { type_search = 3; }
      type_search = 3;
    }    
    if (type_search == 2 &&
	!metaexon->left_closed) {
      type_search = 0;
    } else if (type_search == 3 &&
	       !metaexon->right_closed) {
      type_search = 0;
    }
  } else { 
    //printf("METAEXON NULL!\n");
    type_search = 0; 
  }
  
  
  if (type_search == 2) {
    //Search other
    //printf("RIGHT SINGLE ANCHOR SEARCH\n");
    cigar_code = search_right_single_anchor(read_gap, 
					    cal,
					    0, 
					    metaexon->left_breaks,
					    query,
					    metaexons,
					    genome, 
					    avls_list);    
  } else if (type_search == 3) {
    //printf("Search WITH LEFT ANCHOR\n");
    cigar_code = search_left_single_anchor(read_gap, 
					   cal,
					   0, 
					   metaexon->right_breaks,
					   query,
					   metaexons,
					   genome, 
					   avls_list);   
  } //else {

  
  if (cigar_code == NULL) {
    //printf("Search NORMAL\n");
    if (genome_end - genome_start >= 2048) {
      pthread_mutex_unlock(&metaexons->mutex[cal->chromosome_id - 1]);
      return NULL;
    }    
   
    genome_read_sequence_by_chr_index(reference, 0, 
				      chromosome_id - 1,
				      &genome_start, &genome_end,
				      genome);    
    int distance = 0;
    int err_dsp = 0;
    const int num_err = 2;
    
    //printf("REF: %s\n", reference);
    //printf("FILL EXTREM GAP %i\n", read_gap );
    //printf("Extract [%i:%lu-%lu]: %s\n", chromosome_id, genome_start, genome_end, reference);
    //assert(read_start + read_gap < strlen(query));
    if (read_start + read_gap > strlen(query)) {
      LOG_FATAL_F("%i + %i >= %lu\n", read_start, read_gap , strlen(query));
    }
    //assert(read_gap < strlen(reference));
    for (int i = 0; i < read_gap; i++) {
      //printf("%c vs %c \n", query[read_start + i], reference[i]);
      // i, read_gap, distance);      
      if (query[read_start + i] != reference[i]) {
	distance++;
      }
    }
    
    //printf("distance = %i < max_err = %i\n", distance, ((read_gap / 4) + 1));

    const int MAX_GAP = 10;
    const int MAX_ERR = 3;
   
    assert(distance <= read_gap);
    if ((read_gap <= MAX_GAP && distance <= MAX_ERR) || 
	(read_gap > MAX_GAP && distance <= ((read_gap / 7))) ) {
      cigar_code = cigar_code_new();
      cigar_code->distance = distance;
      cigar_code_append_new_op(read_gap, 'M', cigar_code);
    }
    
  }
  
  pthread_mutex_unlock(&metaexons->mutex[cal->chromosome_id - 1]);


  return cigar_code;

}

int extend_extrem_nt(char *reference, char *query, int extrem_type, int *distance) {
  int pos, action, limit, num_errors = 0;
  int limit_errors = strlen(query)/4, num_matches = 0;
  int consecutive_err = 0;
  int max_consecutive_err = 3;
  int number_nt = 0;

  //printf("Reference: %s\n", reference);
  //printf("Query    : %s\n", query);

  if (extrem_type == LEFT_EXTREM) {
    pos = 0;
    action = 1;
    limit = strlen(query);
  } else {
    pos = strlen(query) - 1;
    action = -1;
    limit = -1;
  }

  while (pos != limit) {
    if (reference[pos] != query[pos]) { 
      num_errors++;
      consecutive_err++;
      if (num_errors >= limit_errors || 
	  consecutive_err > max_consecutive_err) { number_nt = 0; break; }
    } else {
      num_matches++;
      consecutive_err = 0;
    }
    pos += action;
    number_nt++;
  }
  
  //printf("I come to pos %i with %i errors\n", number_nt, num_errors);
  *distance = num_errors;

  return number_nt;

}


float generate_cals_score(array_list_t *cals_list, int read_length) {
  int len_cal = 0;
  size_t num_cals = array_list_size(cals_list);
  cal_t *cal;
  linked_list_iterator_t itr;  
  seed_region_t *s, *s_prev;
  int i;

  //printf("==== CALS SCORE ====\n");
  for (i = 0; i < num_cals; i++) {    
    cal = array_list_get(i, cals_list);
    //printf("\tCAL (%i)[%i:%lu - %lu]:\n", cal->strand, cal->chromosome_id, cal->start, cal->end);
    // s_prev = NULL;
    linked_list_iterator_init(cal->sr_list, &itr);
    s = (seed_region_t *) linked_list_iterator_curr(&itr);
    while (s != NULL) {
      if (s->read_start > s->read_end) { assert(s->read_start); }
      len_cal += s->read_end - s->read_start;
      //if (s_prev) {
      //len_cal += s->read_start - s_prev->read_end;
      //}
      //s_prev = s;
      //printf("\t\tSEED %i:[%lu|%i-%i|%lu]Len CAL %i\n", i, s->genome_start, s->read_start, s->read_end, s->genome_end, len_cal);
      s = (seed_region_t *) linked_list_iterator_next(&itr);
    }

    if (array_list_size(cal->candidates_seeds_start)) {
      len_cal += 16;
      //printf("\tSeeds Candidates %i:\n", array_list_size(cal->candidates_seeds_start));
      for (int j = 0; j < array_list_size(cal->candidates_seeds_start); j++) {
	seed_region_t *seed_region = array_list_get(j, cal->candidates_seeds_start);
	//printf("\t\t\t Candidate Seed S:[%lu|%i-%i|%lu]\n", seed_region->genome_start, seed_region->read_start, seed_region->read_end, seed_region->genome_end);
      }
    }

    if (array_list_size(cal->candidates_seeds_end)) {
      len_cal += 16;
      //printf("\tSeeds Candidates %i:\n", array_list_size(cal->candidates_seeds_end));
      for (int j = 0; j < array_list_size(cal->candidates_seeds_end); j++) {
	seed_region_t *seed_region = array_list_get(j, cal->candidates_seeds_end);
	//printf("\t\t\t Candidate Seed E:[%lu|%i-%i|%lu]\n", seed_region->genome_start, seed_region->read_start, seed_region->read_end, seed_region->genome_end);
      }
    }
    //printf("\t<------- NEW CAL --------\n");
  }
  //printf("(SCORE %i,%i, %f)\n", len_cal, read_length, (float)(len_cal*100)/(float)read_length);
  //printf("<##### END FUNCTION #####>\n");

  return (float)(len_cal*100)/(float)read_length;

}

void order_cals(array_list_t *cals_list) {

  if (array_list_size(cals_list) <= 1) { return; }
 
  cal_t *cal_prev, *cal_next;
  for (int i = 0; i < array_list_size(cals_list) - 1; i++) {
    cal_prev = (cal_t *)array_list_get(i, cals_list);
    for (int j = i + 1; j < array_list_size(cals_list); j++) {
      cal_next = (cal_t *)array_list_get(j, cals_list);
      if (cal_next->strand < cal_prev->strand) {
	//printf("(%i-%i) %i[%i:%lu-%lu] <(strand)> %i[%i:%lu-%lu]\n", i, j, cal_prev->strand, cal_prev->chromosome_id, cal_prev->start, cal_prev->end, 
	//       cal_next->strand, cal_next->chromosome_id, cal_next->start, cal_next->end);
	array_list_swap(i, j, cals_list);
	cal_prev = cal_next;
      } else {
	if (cal_next->chromosome_id < cal_prev->chromosome_id) {
	  //printf("(%i-%i) %i[%i:%lu-%lu] <(chr)> %i[%i:%lu-%lu]\n", i, j, cal_prev->strand, cal_prev->chromosome_id, cal_prev->start, cal_prev->end, 
	  //	 cal_next->strand, cal_next->chromosome_id, cal_next->start, cal_next->end);
	  array_list_swap(i, j, cals_list);
	  cal_prev = cal_next;
	} else {
	  if (cal_next->start < cal_prev->start) {
	    //printf("(%i-%i) %i[%i:%lu-%lu] <(start)> %i[%i:%lu-%lu]\n", i, j, cal_prev->strand, cal_prev->chromosome_id, cal_prev->start, cal_prev->end, 
	    //	   cal_next->strand, cal_next->chromosome_id, cal_next->start, cal_next->end);
	    array_list_swap(i, j, cals_list);
	    cal_prev = cal_next;
	  }
	}
      }
    }
  }

}


int merge_and_filter_cals(array_list_t *cals_targets, 
			  array_list_t *cals_list, 
			  fastq_read_t *fq_read, 
			  bwt_optarg_t *bwt_optarg,
			  bwt_index_t *bwt_index, 
			  genome_t *genome, 
			  float *score_ranking) {

  int num_cals = array_list_size(cals_list);
  int cal_pos;
  cal_t *cal_prev, *cal_next, *cal;
  array_list_t *merge_cals, *fusion_cals;
  seed_region_t *s, *s_prev, *s_next;
  float *cals_score = score_ranking;  
  int number_of_best;
  size_t max_intron_size = 500000;
  float score;
  cigar_code_t *cigar_code = NULL;
  linked_list_iterator_t itr;  
  linked_list_t *linked_list;
  int seed_err_size = 20;
  char query[fq_read->length];
  char reference[fq_read->length];
  char *rev_comp = NULL;
  char *sequence;
  int number_nt;
  size_t genome_start, genome_end;
  int best_cals = num_cals;
  int distance = 0;

  const float MIN_CAL_SCORE = 60.0;
  register int i, j;

  if (!num_cals) { return 0; }

  //===== Step-1: Concatenate CALs =====//
  LOG_DEBUG("STEP-1: CONCATENATE CALS");
  merge_cals = array_list_new(100,
			      1.25f,
			      COLLECTION_MODE_ASYNCHRONIZED);
  cal_pos = 0;
  cal_prev = (cal_t *)array_list_get(cal_pos++, cals_list);
  
  array_list_insert(cal_prev, merge_cals);    
  while (cal_pos < num_cals) {
    cal_next = (cal_t *)array_list_get(cal_pos, cals_list);      
    s_prev = linked_list_get_last(cal_prev->sr_list);
    s = linked_list_get_first(cal_next->sr_list);
    
    assert(s_prev != NULL);
    assert(s != NULL);
    if (cal_prev->chromosome_id == cal_next->chromosome_id && 
	cal_prev->strand == cal_next->strand && 
	s_prev->read_end <= s->read_start &&
	(cal_next->start <= (cal_prev->end + max_intron_size))) {
      //printf("Merge!! cal_prev->end = %lu, cal_next->start = %lu\n", cal_prev->end, cal_next->start);
      array_list_insert(cal_next, merge_cals);
    } else { 
      array_list_insert(merge_cals, cals_targets);
      merge_cals = array_list_new(10,
				  1.25f,
				  COLLECTION_MODE_ASYNCHRONIZED);
      array_list_insert(cal_next, merge_cals);
    }                                                                                         
    cal_prev = cal_next;
    cal_pos++;
  }
  array_list_insert(merge_cals, cals_targets);         
  array_list_clear(cals_list, (void *)NULL);
  //===== Step-1: END =====//


  //===== Step-2: Generate CALs Score =====//
  LOG_DEBUG("STEP-2: GENERATE CALS SCORE");

  if (array_list_size(cals_targets) > 100) {
    LOG_FATAL("MAX CALS OVERFLOW\n");
  }

  for (i = 0; i < array_list_size(cals_targets); i++) {
    fusion_cals = array_list_get(i, cals_targets);
    cals_score[i] = generate_cals_score(fusion_cals, fq_read->length);
    //printf("SCORE %f\n", cals_score[i]);
  }
  
  //===== Step-2: END =====//

  /*num_cals = array_list_size(cals_targets);
  if (best_cals > num_cals) {
    number_of_best = num_cals;
  } else {
    number_of_best = best_cals;
    }*/

  number_of_best = array_list_size(cals_targets);
  
  //===== Step-3: Ranking CALs by score (the n best)=====//
  LOG_DEBUG("STEP-3: ORDER CALS");
  for (i = 0; i < number_of_best; i++) {
    for (j = i + 1; j < array_list_size(cals_targets); j++) {
      if (cals_score[j] > cals_score[i]) {
	array_list_swap(i, j, cals_targets);
	score = cals_score[j];
	cals_score[j] = cals_score[i];
	cals_score[i] = score;
      }
    }
  }
  //===== Step-3: END =====//

  //===== Step-4: Search Near Seeds of CALs =====//
  for (i = 0; i < number_of_best; i++) {
    fusion_cals = array_list_get(i, cals_targets);
    // Has this CAL one or more seeds near?
    cal_t *cal_prev = array_list_get(0, fusion_cals);
    cal_t *cal_next = array_list_get(array_list_size(fusion_cals) - 1, fusion_cals);
    
    //TODO: FIRST SEED? NO... SEARCH THE BEST SEEDS
    if (array_list_size(cal_prev->candidates_seeds_start) > 0) {
      seed_region_t *seed_region = array_list_remove_at(0, cal_prev->candidates_seeds_start);
      assert(seed_region != NULL);
      linked_list_t *linked_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
      linked_list_insert_first(seed_region, linked_list);

      cal = cal_new(cal_prev->chromosome_id, cal_prev->strand, 
		    seed_region->genome_start, seed_region->genome_end, 
		    1, linked_list, 
		    linked_list_new(COLLECTION_MODE_ASYNCHRONIZED));
      array_list_insert_at(0, cal, fusion_cals);
    }

    if (array_list_size(cal_next->candidates_seeds_end) > 0) {
      seed_region_t *seed_region = array_list_remove_at(0, cal_next->candidates_seeds_end);
      assert(seed_region != NULL);
      linked_list_t *linked_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
      linked_list_insert_first(seed_region, linked_list);

      cal = cal_new(cal_next->chromosome_id, cal_next->strand, 
		    seed_region->genome_start, seed_region->genome_end, 
		    1, linked_list,
		    linked_list_new(COLLECTION_MODE_ASYNCHRONIZED));
      array_list_insert(cal, fusion_cals);
    }

    /*for (j = 0; j < array_list_size(fusion_cals); j++) {
      cal = array_list_get(j, fusion_cals);            
      array_list_insert(cal, cals_list);
    }
    */
  }
  
  //Extend extrem gaps to final    
  /*for (i = 0; i < number_of_best; i++) {
    fusion_cals = array_list_get(i, cals_targets);
    cal_prev = array_list_get(0, fusion_cals);
    cal_next = array_list_get(array_list_size(fusion_cals) - 1, fusion_cals);

    s_prev = linked_list_get_first(cal_prev->sr_list);
    s_next = linked_list_get_last(cal_next->sr_list);

    int nt_tot = 0;
    //printf("Extend [%i:%lu-%lu]??\n", cal_prev->chromosome_id, cal_prev->start, cal_next->end);
    if (s_prev->read_start > 20) {
      //printf("\tExtend left extrem\n");
      if (cal_prev->strand == 1) {
	if (!rev_comp ) {
	  rev_comp = (char *) calloc(fq_read->length + 1, sizeof(char));
	  strcpy(rev_comp, fq_read->sequence);
	  seq_reverse_complementary(rev_comp, fq_read->length);
	}
	sequence = rev_comp;
      } else {
	sequence = fq_read->sequence;
      }

      genome_start = cal_prev->start - s_prev->read_start + 1;
      genome_end = cal_prev->start;
      genome_read_sequence_by_chr_index(reference, 0, 
					cal_prev->chromosome_id - 1, &genome_start, &genome_end, genome);      
      memcpy(query, sequence, s_prev->read_start);
      query[s_prev->read_start] = '\0';
      number_nt = extend_extrem_nt(reference, query, LEFT_EXTREM, &distance);
      if (number_nt == 0) {
	q[*num_sw] = strdup(query);
	r[*num_sw] = strdup(reference);
	fusion_coords[(*num_sw)++] = fusion_coords_new(0, 0, 0, 0, 0, 0, 0, 
						       0, FIRST_SW, fq_read->id, cal_prev, read_id, i);
	//printf("Insert Reference CAL %i, Read %i Start\n", i, read_id);
      } else {
	s_prev->read_start -= number_nt;
	s_prev->genome_start -= number_nt;
	cal_prev->start -= number_nt;      
	nt_tot += number_nt;
	cals_score[i] += (((nt_tot - distance)*100)/fq_read->length);
	//printf("Actualization score\n");
      }
    }
    
    if (s_next->read_end < fq_read->length - 20) {
      //printf("\tExtend right extrem\n");
      if (cal_next->strand == 1) {
	if (!rev_comp ) {
	  rev_comp = (char *) calloc(fq_read->length + 1, sizeof(char));
	  strcpy(rev_comp, fq_read->sequence);
	  seq_reverse_complementary(rev_comp, fq_read->length);
	}
	sequence = rev_comp;
      } else {
	sequence = fq_read->sequence;
      }

      genome_start = cal_next->end;
      genome_end = cal_next->end + (fq_read->length - s_next->read_end) - 1;
      genome_read_sequence_by_chr_index(reference, 0, 
					cal_next->chromosome_id - 1, &genome_start, &genome_end, genome);

      //printf("#############From %i:%lu-%lu: %s\n", cal_next->chromosome_id - 1, genome_start, genome_end, reference);

      memcpy(query, sequence + s_next->read_end, fq_read->length - s_next->read_end);
      query[fq_read->length - s_next->read_end] = '\0';
      number_nt = extend_extrem_nt(reference, query, RIGHT_EXTREM, &distance);
      if (number_nt == 0) {
	q[*num_sw] = strdup(query);
	r[*num_sw] = strdup(reference);
	fusion_coords[(*num_sw)++] = fusion_coords_new(0, 0, 0, 0, 0, 0, 0, 0, 
						       LAST_SW, fq_read->id, cal_next, read_id, i);
	//printf("Insert Reference CAL %i, Read %i End\n", i, read_id);
      } else {
	s_next->read_end += number_nt;
	s_next->genome_end += number_nt;
	cal_next->end += number_nt;	
	nt_tot += number_nt;      
	//printf("Actualization score %i nt, %i errors\n", number_nt, distance);
	cals_score[i] += (((nt_tot - distance)*100)/fq_read->length);
      }
    }    
  }*/
  
  //The best CAL has a start/end big gap?
  //printf("BEST SCORE %f\n", cals_score[0]);
  /*
    if (cals_score[0] <= MIN_CAL_SCORE) {
    fusion_cals = array_list_get(0, cals_targets);
    cal = array_list_get(0, fusion_cals);
    s_prev = linked_list_get_first(cal->sr_list);
    s_next = linked_list_get_last(cal->sr_list);
  
    //printf("$$$Make Seeds with one Error!! FIRST SEED:[%i-%i] LAST_SEED:[%i-%i] %s\n", s_prev->read_start, s_prev->read_end,
    //	   s_next->read_start, s_next->read_end, fq_read->id);
    
    if (s_prev->read_start >= seed_err_size) {
      printf("\t @@@@SEEDs in First positions [%i-%i]\n", 0, s_prev->read_start - 1);
      //Seeds in Start position
      array_list_t *mapping_list_prev = array_list_new(1000,
						       1.25f,
						       COLLECTION_MODE_ASYNCHRONIZED);
      
      bwt_map_inexact_seeds_by_region(0, s_prev->read_start - 1,
				      cal->chromosome_id, cal->start - 500000,
				      cal->start,
				      fq_read->sequence, seed_err_size,
				      seed_err_size,
				      bwt_optarg,
				      bwt_index,
				      mapping_list_prev);
      
      
      for (int r = 0; r < array_list_size(mapping_list_prev); r++) {
	region_t *region = array_list_get(r, mapping_list_prev);
	printf("\t Region [%i:%lu-%lu]\n", region->chromosome_id, region->start, region->end);
      }

      array_list_free(mapping_list_prev, region_bwt_free);

    }

    if (s_next->read_end <= fq_read->length - seed_err_size) {
      printf("\t @@@@SEEDs in Last positions [%i-%i]\n", s_next->read_end, fq_read->length - 1);
      //Seeds in End position
      array_list_t *mapping_list_next = array_list_new(1000,
						       1.25f,
						       COLLECTION_MODE_ASYNCHRONIZED);
     
      bwt_map_inexact_seeds_by_region(s_next->read_end, fq_read->length - 1,
				      cal->chromosome_id, cal->end,
				      cal->end + 500000,
				      fq_read->sequence, seed_err_size,
				      seed_err_size/2,
				      bwt_optarg,
				      bwt_index,
				      mapping_list_next);
      
      for (int r = 0; r < array_list_size(mapping_list_next); r++) {
	region_t *region = array_list_get(r, mapping_list_next);
	printf("\t Region [%i:%lu-%lu]\n", region->chromosome_id, region->start, region->end);
      }

      array_list_free(mapping_list_next, region_bwt_free);
    }
    }*/
  //End

  //Delete other CALs
  /*for (i = array_list_size(cals_targets) - 1; i >= number_of_best; i--) {
    fusion_cals = array_list_remove_at(i, cals_targets);
    for (j = 0; j < array_list_size(fusion_cals); j++) {
      cal = array_list_get(j, fusion_cals);
      cigar_code = (cigar_code_t *)cal->info;
      if (cigar_code != NULL) {
	array_list_clear(cigar_code->ops, cigar_op_free);
	cigar_code_free(cigar_code);
      }
      cal_free(cal);
    }
    array_list_free(fusion_cals, NULL);
  }
  */
  //free(rev_comp);

  return number_of_best;

}




void generate_reference_splice_juntion(array_list_t *cals_targets, char *query_revcomp, 
				       fastq_read_t *fq_read, char **r, char **q, 
				       fusion_coords_t **fusion_coords, int *num_sw, 
				       genome_t *genome, size_t id_read) {
  int n_fusion_cals = array_list_size(cals_targets);

  cal_t *cal_prev, *cal_next, *cal;
  array_list_t *fusion_cals;
  char *query_map;
  seed_region_t *s, *s_prev;
  char reference[2048];
  char reference_prev[2048];
  char reference_next[2048];
  char query[2048];

  size_t genome_start, genome_end;
  int read_start, read_end;
  size_t genome_start2, genome_end2;
  int seeds_nt;
  int flank = 30;

  register int i, j;

  //Delete for debuging, detect splice junctions
  for (i = 0; i < n_fusion_cals; i++) {
    fusion_cals = array_list_get(i, cals_targets);
    j = 0;
    cal_prev = array_list_get(j, fusion_cals);
    s_prev = linked_list_get_last(cal_prev->sr_list);
    if (cal_prev->strand == 1) {
      query_map = query_revcomp;
    } else {
      query_map = fq_read->sequence;
    }

    //printf("SP_CAL_PREV(%i) %i:%lu-%lu %i:%i=%i\n", cal->strand, cal->chromosome_id, cal_prev->end, 
    //     cal->start, read_start, read_end, seeds_nt);
    
    for (j = 1; j < array_list_size(fusion_cals); j++) {
      cal = array_list_get(j, fusion_cals);
      s = linked_list_get_first(cal->sr_list);

      read_start = s_prev->read_end;
      read_end = s->read_start;
      seeds_nt = read_end - read_start;
      
      //printf("SP_CAL(%i) %i:%lu-%lu %i:%i=%i\n", cal->strand, cal->chromosome_id, cal_prev->end, 
      //   cal->start, read_start, read_end, seeds_nt);
      
      if (seeds_nt >= 10) {
	//Extend to Right --> <-- Extend to Left
	genome_start = s_prev->genome_end;
	genome_end = s_prev->genome_end + seeds_nt - 1;
	genome_read_sequence_by_chr_index(reference_prev, 0, 
					  cal->chromosome_id - 1, &genome_start, &genome_end, genome);

	genome_start2 = s->genome_start - seeds_nt;
	genome_end2 = s->genome_start - 1;
	genome_read_sequence_by_chr_index(reference_next, 0, 
					  cal->chromosome_id - 1, &genome_start2, &genome_end2, genome);

	memcpy(query, query_map + read_start, read_end - read_start);
	query[read_end - read_start]  = '\0';

	//CAll new function
	int dsp_e1, dsp_e2;
	int lim_err = 3;

	extend_by_mismatches(reference_prev, reference_next, query, 
			     0, strlen(reference_next) - 1, 
			     0, read_end - read_start - 1, lim_err,
			     &dsp_e1, &dsp_e2);
	//printf("dsp_e1=%i, dsp_e2=%i\n", dsp_e1, dsp_e2);

	cal_prev->end = cal_prev->end + dsp_e1;	
	cal->start = cal->start - dsp_e2;

	s_prev->read_end = s_prev->read_end + dsp_e1;
	s->read_start = s->read_start - dsp_e2;

	if (s_prev->read_end > s->read_start) { 
	  seeds_nt = s_prev->read_end - s->read_start;
	  cal_prev->end -= seeds_nt;
	  cal->start += seeds_nt;
	  s_prev->read_end-= seeds_nt;
	  s->read_start += seeds_nt;
	}
      }
      
      read_start = s_prev->read_end - flank;
      if (read_start < 0) { read_start = 0; }
      read_end = s->read_start + flank;
      if (read_end >= fq_read->length) { read_end = fq_read->length - 1; }

      //Extract and fusion Reference SW
      genome_start = cal_prev->end - flank;
      genome_end = cal_prev->end + flank - 1;
      genome_read_sequence_by_chr_index(reference_prev, 0, 
					cal_prev->chromosome_id - 1, &genome_start, &genome_end, genome);

      //printf("Ref Prev %i\n", strlen(reference_prev));
      genome_start2 = cal->start - flank;
      genome_end2 = cal->start + flank - 1;
      genome_read_sequence_by_chr_index(reference_next, 0, 
					cal->chromosome_id - 1, &genome_start2, &genome_end2, genome);
      
      //printf("Ref Next %i [%i:%i-%i]=%s\n", strlen(reference_next), cal->chromosome_id, genome_start2, genome_end2, reference_next);
      //if (genome_start2 < genome_start) { printf("%s\n", fq_read->id); assert(genome_start); }
      //printf("\tSP_CAL(%i) %i:%lu-%lu %i:%i\n", cal->strand, cal->chromosome_id, cal_prev->end, 
      //   cal->start, read_start, read_end);
      
      strcat(reference_prev, reference_next);

      if (read_start > read_end) { printf("ALERT ERROR:: %s\n", fq_read->id);exit(-1); }
      //printf("Read %i-%i: %s\n", read_start, read_end, query_map);

      memcpy(query, query_map + read_start, read_end - read_start);
      query[read_end -read_start] = '\0';

      r[*num_sw] = strdup(reference_prev);
      q[*num_sw] = strdup(query);

      //printf("Que: %s\n", q[*num_sw]);
      //printf("Ref: %s\n", r[*num_sw]);
      //printf("Ref len %i, Query len %i\n", strlen(r[num_sw]), strlen(q[num_sw]));
      fusion_coords[(*num_sw)++] = fusion_coords_new(genome_start, genome_end,
						     genome_start2, genome_end2, 
						     read_start, read_end,
						     cal->chromosome_id, cal->strand,
						     MIDDLE_SW, fq_read->id, cal_prev, id_read, 0);	
      
      cal_prev = cal;
      s_prev = linked_list_get_last(cal_prev->sr_list);
      
      }
  }
}


/*
array_list_t *fusion_regions (array_list_t *regions_list, int max_distance) {
  int r = 0;
  region_t *region, *region_next;
  array_list_t *merge_regions_list = array_list_new(array_list_size(regions_list),
						   1.25f,
						   COLLECTION_MODE_ASYNCHRONIZED);
  int last_id;

  //printf("Regions Found %i:\n", array_list_size(regions_list));
  while (array_list_size(regions_list) > 0) {
    region_t *region = array_list_get(r, regions_list);
    //printf("\t::: %i-Region[%lu|%i-%i|%lu] (%i/%i)\n", 
    //	   region->id, region->start, region->seq_start, region->seq_end, region->end, 
    //	   r, array_list_size(regions_list));
    if (r + 1 < array_list_size(regions_list)) {
      region_t *region_next = array_list_get(r + 1, regions_list);
      if (region->id == region_next->id) {
	//TODO: Equal seeds id select not do this
	r++;
	continue;
      } else if ((region->end + max_distance) >= region_next->start) {
	//Fusion seeds
	region_next->start = region->start;
	region_next->seq_start = region->seq_start;
      } else {
	//Not fusion
	//printf("\t\tInsert region\n");
	array_list_insert(region_bwt_new(region->chromosome_id,
					 region->strand,
					 region->start,
					 region->end,
					 region->seq_start,
					 region->seq_end,
					 region->seq_len,
					 region->id), 
			  merge_regions_list);
      }
    } else {
      //printf("\t\tInsert region\n");
      array_list_insert(region_bwt_new(region->chromosome_id,
				       region->strand,
				       region->start,
				       region->end,
				       region->seq_start,
				       region->seq_end,
				       region->seq_len,
				       region->id),
			merge_regions_list);	    
      break;
    }
    r++;
  }

  array_list_clear(regions_list, region_bwt_free);

  for (int i = 0; i < array_list_size(merge_regions_list); i++) {
    region_t *region = array_list_get(i, merge_regions_list);
    array_list_insert(region, regions_list);
  }
  
  array_list_free(merge_regions_list, NULL);

  return regions_list;
  
}
*/

//============================ STRUCTURES AND TYPES SECTION =============================//



typedef struct sw_item {
  int type_sw;
  int read_id;
  int fusion_id;
  int cal_id;
  seed_region_t *seed_prev;
  seed_region_t *seed_next;
  cal_t *cal_prev;
  cal_t *cal_next;
  meta_alignment_t *meta_alignment;
  void *info;
} sw_item_t;

sw_item_t *sw_item_new(int type_sw, int read_id, 
		       int fusion_id, int cal_id,
		       cal_t *cal_prev, cal_t *cal_next, 
		       meta_alignment_t *meta_alignment, 
		       seed_region_t *seed_prev,
		       seed_region_t *seed_next,
		       void *info);

typedef struct sw_depth {
  char *q[MAX_DEPTH];
  char *r[MAX_DEPTH];
  sw_item_t *items[MAX_DEPTH];
  int depth;
} sw_depth_t;

void sw_depth_insert(char *query, char *reference, sw_item_t *sw_item, 
		     sw_optarg_t *sw_optarg, sw_multi_output_t *output, 
		     avls_list_t *avls_list, metaexons_t *metaexons, 
		     sw_depth_t *sw_depth, genome_t *genome, int min_intron_size);

//=======================================================================================//


//=============================== META ALIGNMENT SECTION ================================//
meta_alignment_t *meta_alignment_new() {
  meta_alignment_t *meta_alignment = (meta_alignment_t *)malloc(sizeof(meta_alignment_t));

  meta_alignment->cals_list = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  meta_alignment->status = META_OPEN;
  meta_alignment->sp_sw = 0;
  meta_alignment->num_cigars = 0;
  meta_alignment->cigar_code = cigar_code_new();
  meta_alignment->type = META_ALIGNMENT_NONE;
  meta_alignment->cigar_left = NULL;
  meta_alignment->cigar_right = NULL;
  meta_alignment->score = 0;
  meta_alignment->flag = 0;

  return meta_alignment;

}

void meta_alignment_free(meta_alignment_t *meta_alignment) {
  //if (meta_alignment->cals_list != NULL) { array_list_free(meta_alignment->cals_list, NULL); }
  cigar_code_free(meta_alignment->cigar_code);
  free(meta_alignment);
}

int meta_alignment_num_cals(meta_alignment_t *meta_alignment) {
  return array_list_size(meta_alignment->cals_list);
}

cal_t *meta_alignment_get_first_cal(meta_alignment_t *meta_alignment) {
  return array_list_get(0, meta_alignment->cals_list);
}

cal_t *meta_alignment_get_last_cal(meta_alignment_t *meta_alignment) {
  return array_list_get(meta_alignment_num_cals(meta_alignment) - 1, meta_alignment->cals_list);
}


int meta_alignment_num_cigars(meta_alignment_t *meta_alignment) {
  if (meta_alignment != NULL) {
    return meta_alignment->num_cigars;
  } else {
    return 0;
  }
}

int meta_alignment_set_status(int status, meta_alignment_t *meta_alignment) {
  meta_alignment->status = status;
  return status;
}

int meta_alignment_get_status(meta_alignment_t *meta_alignment) {
  return meta_alignment->status;
}

int meta_alignment_insert_cal(cal_t *cal, meta_alignment_t *meta_alignment) {
  return array_list_insert(cal, meta_alignment->cals_list);
}


void meta_alignments_order_by_score(array_list_t *meta_alignments) {
  meta_alignment_t *meta_prev, *meta_next;
  for (int i = 0; i < array_list_size(meta_alignments) - 1; i++) {
    meta_prev = array_list_get(i, meta_alignments);
    for (int j = i + 1; j < array_list_size(meta_alignments); j++) {
      meta_next = array_list_get(j, meta_alignments);
      if (meta_next->score > meta_prev->score) {
	array_list_swap(i, j, meta_alignments);
      }
    }
  }
  
}

void meta_alignment_insert_cigar(cigar_code_t *cigar, int type, int pos, meta_alignment_t *meta_alignment) {  
  //if (cigar != NULL && (type == CIGAR_ANCHOR_LEFT || type == CIGAR_ANCHOR_RIGHT)) {
    //printf("-----------------------------> INSERT %s. (%i)\n", new_cigar_code_string(cigar), array_list_size(cigar->ops));
  //}

  if (meta_alignment == NULL) { return; }
  else {
    if (type == CIGAR_ANCHOR_LEFT) {
      meta_alignment->cigar_right = cigar;
    } else if (type == CIGAR_ANCHOR_RIGHT) {
      meta_alignment->cigar_left = cigar;
    } else {
      //printf("Insert middle cigar %i => %s\n", pos, new_cigar_code_string(cigar));
      //if (pos > 20 || pos < 0) { exit(-1); }
      meta_alignment->middle_cigars[pos] = cigar;
      meta_alignment->type_cigars[pos] = type;
      meta_alignment->num_cigars++;
    }
  }
}

void meta_alignment_calculate_score(meta_alignment_t *meta_alignment) {
  cigar_code_t *cigar_code = meta_alignment->cigar_code;
  int num_di = 0;
  meta_alignment->score = 0;
  if (cigar_code == NULL) { printf("NULL CIGAR\n"); return; }
  //printf("NUM OPS: %i\n", array_list_size(cigar_code->ops));
  for (int i = 0; i < array_list_size(cigar_code->ops); i++) {
    cigar_op_t *op = array_list_get(i, cigar_code->ops);
    if (op->name == 'M' || op->name == 'I') { 
      meta_alignment->score += op->number; 
    }
    //if (op->name == 'I' || op->name == 'D') {
    //num_di += op->number;
    //}
    //printf("SCORE (%i) OP+ : %i%c\n", meta_alignment->score, 
    //	   op->number, op->name);
  }
}

void meta_alignment_calculate_f_score(meta_alignment_t *meta_alignment) {
  cigar_code_t *cigar_code = meta_alignment->cigar_code;
  if (cigar_code == NULL) { meta_alignment->f_score = 0; }
  
  
}

int meta_alignment_get_cals_score(meta_alignment_t *meta_alignment) {
  int num_di = 0, score = 0;
  printf("GET SCORE\n");
  for (int i = 0; i < array_list_size(meta_alignment->cals_list); i++) {
    cal_t *cal = array_list_get(i, meta_alignment->cals_list);
    cal_print(cal);
    linked_list_item_t *item = cal->sr_list->first;
    while (item != NULL) {
      seed_region_t *seed_aux = item->item;
      score += seed_aux->read_end - seed_aux->read_start;
      item= item->next;
    }
  }

  return score;

}

void meta_alignment_close(meta_alignment_t *meta_alignment) {
  cal_t *first_cal, *cal;
  cigar_code_t *cigar_code, *cigar_code_aux;
  linked_list_item_t *list_item;
  seed_region_t *s, *s_prev;
  int type;
  int cr_pos = 0;
  cigar_op_t *op;
  int bad_cigar = 0;

  //printf("\n==================== CLOSE META ALIGNMENT ==========================\n");
  assert(meta_alignment_num_cals(meta_alignment) > 0);
  cigar_code = meta_alignment->cigar_code;
  assert(cigar_code != NULL);

  if (cigar_code_get_num_ops(cigar_code) > 0) {
    array_list_clear(cigar_code->ops, (void *)cigar_op_free);
    cigar_code->distance = 0;
  }

  //cal_t *cal_prev = NULL;
  if (meta_alignment->type == META_ALIGNMENT_MIDDLE) {
    //printf("META_ALIGNMENT_MIDDLE %i\n", meta_alignment->type);
    if (meta_alignment->cigar_left != NULL) {
      cigar_code_aux = meta_alignment->cigar_left;
      for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
	cigar_op_t *op = array_list_get(i, cigar_code_aux->ops);
	//printf("\t OP-SP-LEFT: %i%c\n", op->number, op->name);
	cigar_code_append_new_op(op->number, op->name, cigar_code);
      } 
      cigar_code->distance += cigar_code_aux->distance;
    }

    for (int i = 0; i < meta_alignment_num_cals(meta_alignment); i++) {
      cal_t *cal = array_list_get(i, meta_alignment->cals_list);    
      /*printf(" ***** CLOSE GAP *****\n");
      cal_print(cal);
      printf(" *********************\n");*/
      cigar_code_aux = cal->info;

      if (cr_pos > 0 &&
	  meta_alignment->type_cigars[cr_pos - 1] == CIGAR_SW_MIDDLE) {
	cigar_code_t *c_c = cigar_code_new();
	for (int t = 0; t < array_list_size(cigar_code_aux->ops); t++) {
	  op = array_list_get(t, cigar_code_aux->ops);
	  cigar_code_append_new_op(op->number, op->name, c_c);
	}
	cigar_code_delete_nt(cal->r_flank, 0, c_c);	  
	for (int j = 0; j < cigar_code_get_num_ops(c_c); j++) {
	  op = array_list_get(j, c_c->ops);
	  //printf("\t OP-->1: %i%c\n", op->number, op->name);
	  cigar_code_append_op(op, cigar_code);
	}
	cigar_code_free(c_c);
      } else {      
	for (int j = 0; j < cigar_code_get_num_ops(cigar_code_aux); j++) {
	  op = array_list_get(j, cigar_code_aux->ops);
	  //printf("\t OP-->2: %i%c\n", op->number, op->name);
	  cigar_code_append_new_op(op->number, op->name, cigar_code);
	  //cigar_code_append_op(op, cigar_code);
	}
      }

      cigar_code->distance += cigar_code_aux->distance;

      if (cr_pos < meta_alignment_num_cigars(meta_alignment)) {
	if (meta_alignment->type_cigars[cr_pos] == CIGAR_SW_MIDDLE) {
	  cigar_code_delete_nt(cal->l_flank, 1, cigar_code);	  
	}

	//printf("CLOSE-META: MIDDLE SPLICE\n");
	cigar_code_aux = meta_alignment->middle_cigars[cr_pos++];
	if (cigar_code_aux == NULL) { /*printf("\txxxx exit with bad.\n");*/ bad_cigar = 1; break; }

	for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
	  cigar_op_t *op = array_list_get(i, cigar_code_aux->ops);
	  //cigar_code_append_op(op, cigar_code);
	  cigar_code_append_new_op(op->number, op->name, cigar_code);
	  //printf("\t OP-SP-M : %i%c\n", op->number, op->name);
	}
      }
    }

    if (meta_alignment->cigar_right != NULL) {
      cigar_code_aux = meta_alignment->cigar_right;
      for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
	cigar_op_t *op = array_list_get(i, cigar_code_aux->ops);
	//printf("\t OP-SP-RIGHT: %i%c\n", op->number, op->name);
	cigar_code_append_new_op(op->number, op->name, cigar_code);
      }
      cigar_code->distance += cigar_code_aux->distance;
    }
  } else {    
    //printf("META_ALIGNMENT_OTHER %i\n", meta_alignment->type);
    if (meta_alignment->cigar_left != NULL) {
      cigar_code_aux = meta_alignment->cigar_left;
      for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
	cigar_op_t *op = array_list_get(i, cigar_code_aux->ops);
	//printf("\t OP-SP-RIGHT: %i%c\n", op->number, op->name);
	//cigar_code_append_op(op, cigar_code);
	cigar_code_append_new_op(op->number, op->name, cigar_code);
      } 
      cigar_code->distance += cigar_code_aux->distance;
    }

    cal_t *cal = array_list_get(0, meta_alignment->cals_list);    
    //printf(" ***** CLOSE GAP *****\n");
    //cal_print(cal);
    //printf(" *********************\n");
    cigar_code_aux = cal->info;  
    //printf("SEEDS OPS: %i\n", cigar_code_get_num_ops(cigar_code_aux));
    for (int j = 0; j < cigar_code_get_num_ops(cigar_code_aux); j++) {
      op = array_list_get(j, cigar_code_aux->ops);
      //printf("\t OP: %i%c\n", op->number, op->name);
      //cigar_code_append_op(op, cigar_code);
      cigar_code_append_new_op(op->number, op->name, cigar_code);
    } 

    cigar_code->distance += cigar_code_aux->distance;

    if (meta_alignment->cigar_right != NULL) {
      cigar_code_aux = meta_alignment->cigar_right;
      for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
	cigar_op_t *op = array_list_get(i, cigar_code_aux->ops);
	//printf("\t OP-SP-LEFT: %i%c\n", op->number, op->name);
	//cigar_code_append_op(op, cigar_code);
	cigar_code_append_new_op(op->number, op->name, cigar_code);
      } 
      cigar_code->distance += cigar_code_aux->distance;
    } 
  }

  if (meta_alignment->cigar_code == NULL) { bad_cigar = 1; }
  
  if (!bad_cigar) {
    meta_alignment_set_status(META_CLOSE, meta_alignment);
    meta_alignment_calculate_score(meta_alignment);
  } else {
    array_list_clear(cigar_code->ops, (void *)cigar_op_free);
  }
  //printf("----- META CLOSE INSERT %s.\n", new_cigar_code_string(meta_alignment->cigar_code));
  //printf("------------------------------------------------------------------\n");
}

void merge_seeds_cal(cal_t *cal) {
  seed_region_t *seed_prev = NULL, *seed_next;
  linked_list_item_t *list_item = cal->sr_list->first, *list_item_prev;
  //if (list_item == NULL) { 
  //printf("MERGE SEEDS CAL\n");
  //}
  cigar_code_t *cigar_code = cigar_code_new();
  cigar_op_t *op;
  
  while (list_item != NULL) {
    seed_next = (seed_region_t *)list_item->item;
    //printf(" Seed ---> %i-%i\n", seed_next->read_start, seed_next->read_end);
    if (seed_prev != NULL) {
      if (seed_prev->fusion_right == 1 &&
	  seed_next->fusion_left == 1) {

	if (seed_prev->info == NULL) {
	  //Detect splice
	  //printf("1.Merge %i\n", cal->type_seeds);
	  /*if (cal->type_seeds && array_list_size(seed_prev->errors_list)) {
	    //printf("1.Merge CIGARs....\n");
	    int num_err = array_list_size(seed_prev->errors_list);
	    //Order Errsors cigar
	    for (int k = 0; k < num_err; k++) {
	      bwt_err_t *err_prev = array_list_get(k, seed_prev->errors_list);
	      for (int k1 = k + 1; k1 < num_err; k1++) {
		bwt_err_t *err_next = array_list_get(k1, seed_prev->errors_list);
		if (err_prev->pos > err_next->pos) {
		  array_list_swap(k, k1, seed_prev->errors_list);
		}
	      }
	    }
	    
	    bwt_err_t *err_prev;
	    //Create cigar
	    for (int k = 0; k < num_err; k++) {
	      err_prev = array_list_get(k, seed_prev->errors_list);
	      op = cigar_op_new(err_prev->pos - seed_prev->read_start, 'M');
	      cigar_code_append_op(op, cigar_code);	    
	      
	      cigar_code_append_new_op(1, err_prev->name, cigar_code);
	    }
	    op = cigar_op_new(seed_prev->read_end - err_prev->pos + 1, 'M');
	    cigar_code_append_op(op, cigar_code);	    
	    
	    } else {*/
	    //printf("Noooooooooo 1\n");
	  op = cigar_op_new(seed_prev->read_end - seed_prev->read_start + 1, 'M');
	  cigar_code_append_op(op, cigar_code);
	  //}
	} else {
	  //printf("1.Merge %i\n", cal->type_seeds);
	  cigar_code_t *cigar_code_aux = seed_prev->info;
	  for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
	    op = array_list_get(i, cigar_code_aux->ops);
	    cigar_code_append_op(op, cigar_code);
	  }
	  cigar_code->distance += cigar_code_aux->distance;
	  cigar_code_free(cigar_code_aux);
	  seed_prev->info = NULL;
	}
      } else {
	//Detect splice
	//printf("2.Merge %i\n", cal->type_seeds);
	/*if (cal->type_seeds && array_list_size(seed_prev->errors_list)) {
	  //printf("2.Merge CIGARs....\n");
	  int num_err = array_list_size(seed_prev->errors_list);
	  //Order Errsors cigar
	  for (int k = 0; k < num_err; k++) {
	    bwt_err_t *err_prev = array_list_get(k, seed_prev->errors_list);
	    for (int k1 = k + 1; k1 < num_err; k1++) {
	      bwt_err_t *err_next = array_list_get(k1, seed_prev->errors_list);
	      if (err_prev->pos > err_next->pos) {
		array_list_swap(k, k1, seed_prev->errors_list);
	      }
	    }
	  }
	  
	  bwt_err_t *err_prev;
	  //Create cigar
	  for (int k = 0; k < num_err; k++) {
	    err_prev = array_list_get(k, seed_prev->errors_list);
	    op = cigar_op_new(err_prev->pos - seed_prev->read_start, 'M');
	    cigar_code_append_op(op, cigar_code);	    
	    
	    cigar_code_append_new_op(1, err_prev->name, cigar_code);
	  }
	  op = cigar_op_new(seed_prev->read_end - err_prev->pos + 1, 'M');
	  cigar_code_append_op(op, cigar_code);	    

	  } else {*/
	  //printf("2---.Merge %i\n", cal->type_seeds);
	op = cigar_op_new(seed_prev->read_end - seed_prev->read_start + 1, 'M');
	cigar_code_append_op(op, cigar_code);
	//}
	//printf("MINI SPLICE %lu - %lu = %lu\n", seed_prev->genome_end, 
	//     seed_next->genome_start, 
	//     seed_next->genome_start - seed_prev->genome_end + 1);
	op = cigar_op_new(seed_next->genome_start - seed_prev->genome_end - 1, 'N');
	cigar_code_append_op(op, cigar_code);
      } 
    }
    seed_prev = seed_next;
    list_item = list_item->next;
  }

  if (seed_prev == NULL) { printf("seed prev NULL\n"); exit(-1); }

  /*if (cal->type_seeds && array_list_size(seed_prev->errors_list)) {
    //printf("3.Merge CIGARs....\n");
    int num_err = array_list_size(seed_prev->errors_list);
    //Order Errsors cigar
    for (int k = 0; k < num_err; k++) {
      bwt_err_t *err_prev = array_list_get(k, seed_prev->errors_list);
      for (int k1 = k + 1; k1 < num_err; k1++) {
	bwt_err_t *err_next = array_list_get(k1, seed_prev->errors_list);
	if (err_prev->pos > err_next->pos) {
	  array_list_swap(k, k1, seed_prev->errors_list);
	}
      }
    }
    bwt_err_t *err_prev;
    int start = seed_prev->read_start;
    //Create cigar
    for (int k = 0; k < num_err; k++) {
      err_prev = array_list_get(k, seed_prev->errors_list);
      op = cigar_op_new(err_prev->pos - start, 'M');
      cigar_code_append_op(op, cigar_code);	    
      
      //printf("Err : %c, %i (Add op %i%c1%c)\n", err_prev->name, err_prev->pos, op->number, op->name, err_prev->name);
      cigar_code_append_new_op(1, err_prev->name, cigar_code);
      start = err_prev->pos + 1;
    }

    //printf("Err : %iM\n", seed_prev->read_end - err_prev->pos + 1);
    op = cigar_op_new(seed_prev->read_end - err_prev->pos + 1, 'M');
    cigar_code_append_op(op, cigar_code);
    } else {    */
    //printf("3.Merge CIGARs --- \n");
  op = cigar_op_new(seed_prev->read_end - seed_prev->read_start + 1, 'M');
  cigar_code_append_op(op, cigar_code);
  //}
  
  //printf("FILL GAPS CLOSE CAL [%i:%lu-%lu]: %s\n", cal->chromosome_id, cal->start, cal->end, new_cigar_code_string(cigar_code));

  cal->info = cigar_code;
  
}

void meta_alignment_fill_gaps(int meta_type,
			      meta_alignment_t *meta_alignment, 
			      char *query_map,
			      genome_t *genome,
			      sw_optarg_t *sw_optarg,
			      sw_multi_output_t *output,			     
			      metaexons_t *metaexons, 
			      sw_depth_t *sw_depth,
			      avls_list_t *avls_list, 
			      int min_intron_size) {
  cigar_code_t *cigar_code;
  int max_size = 2048;
  char *reference = (char *)calloc(max_size, sizeof(char));
  int distance;
  int first, last;
  int min_distance;
  int closed;  
  avl_node_t *node_avl_start;
  avl_node_t *node_avl_end;
  int add_seed, sw_add, sp_found;
  sw_item_t *sw_item;
  const int FLANK = 2;	    

  meta_alignment->type = meta_type;
  cigar_code = meta_alignment->cigar_code;

  //printf("FILL GAPS TOTAL CALS %i\n", array_list_size(meta_alignment->cals_list));
  for (int i = 0; i < array_list_size(meta_alignment->cals_list); i++) {
    //printf("Process CAL %i\n", i);
    cal_t *cal = array_list_get(i, meta_alignment->cals_list);
    //printf("===== FILL GAPS CAL: =====\n");
    //cal_print(cal);
    //printf("=========================\n");
    //fill gaps
    cal->num_targets = 0;
    seed_region_t *seed_prev = NULL, *seed_next;
    linked_list_item_t *list_item = cal->sr_list->first, *list_item_prev;    
    while (list_item != NULL) {
      seed_next = (seed_region_t *)list_item->item;
      int num_match = seed_next->read_end - seed_next->read_start + 1;
      if (seed_prev != NULL) {
	//Close nt
	assert(seed_next->read_start > seed_prev->read_end);
	assert(seed_next->genome_start > seed_prev->genome_end);
	size_t gap_read = seed_next->read_start - seed_prev->read_end  - 1;
	//assert(gap_read >= 0);
	size_t gap_genome = seed_next->genome_start - seed_prev->genome_end  - 1; 
	//assert(gap_genome > 0);

	//printf("FILL GAPS ::: (gap_read = %lu), (gap_genome = %lu)\n", gap_read, gap_genome);
	distance = 0;
	closed = 0;
	add_seed = 0;
	sw_add = 0;
	sp_found = 0;
	if (gap_read == gap_genome &&
	    gap_read == 0) {
	  closed = 1;
	} else if (gap_read == gap_genome) {
	  size_t genome_end   = seed_next->genome_start - 1;
          size_t genome_start = seed_prev->genome_end + 1;
	  int read_end        = seed_next->read_start - 1;
	  int read_start      = seed_prev->read_end + 1;
	  char *query         = &query_map[read_start];
	  first = -1;
	  last = -1;
	  //assert(genome_end - genome_start + 1 < 2048);
	  if (genome_end - genome_start >= max_size) {
	    free(reference);
	    max_size = genome_end - genome_start + 1024;
	    reference = (char *)calloc(max_size, sizeof(char));
	  }
          genome_read_sequence_by_chr_index(reference, 0, cal->chromosome_id - 1,
                                            &genome_start, &genome_end, genome);
	  //printf("[%lu|%i]-GAP-[%i|%lu]: %s\n", genome_start, read_start, read_end, genome_end, reference);
	  for (int k = 0; k < gap_read; k++) {
	    //printf("[q:%c vs r:%c]\n", query[k], reference[k]);
	    if (query[k] != reference[k]) {
	      distance++;
	      if (first == -1) first = k;
	      last = k;
	    }
	  }
	  min_distance = (gap_genome / 3) + 2;
	  //printf("Distance %i <= %i\n", distance, min_distance);
	  if (distance <= min_distance) {
	    cigar_code->distance += distance;
	    closed = 1;
	    add_seed = 1;
	  }
	}

	if (!closed) {
	  int gap = gap_genome - gap_read; 
	  if (gap > 40) {
	    //printf("Search splice");
	    int distance_aux;
	    size_t sp_start, sp_end;
	    int sp_type;
	    int nt = search_simple_splice_junction(seed_prev, seed_next,
						   cal->chromosome_id, cal->strand, 
						   query_map, genome, 
						   &sp_start, &sp_end,
						   &sp_type,
						   &distance_aux);

	    if (nt) {
	      cigar_code->distance = distance_aux;
	      int sp_strand = (sp_type == CT_AC_SPLICE ? 1 : 0);

	      allocate_start_node(cal->chromosome_id - 1,
				  sp_strand,
				  sp_start,
				  sp_end,
				  sp_start,
				  sp_end,
				  FROM_READ,
				  sp_type,
				  NULL, 
				  &node_avl_start,
				  &node_avl_end, 
				  avls_list);		  

	      assert(seed_prev->genome_start < node_avl_start->position);
	      assert(node_avl_end->position < seed_next->genome_end);
	      
	      metaexon_insert(cal->strand, cal->chromosome_id - 1,
			      seed_prev->genome_start, node_avl_start->position, 40,
			      METAEXON_RIGHT_END, node_avl_start,
			      metaexons);
	      
	      metaexon_insert(cal->strand, cal->chromosome_id - 1,
			      node_avl_end->position, seed_next->genome_end, 40,
			      METAEXON_LEFT_END, node_avl_end,
			      metaexons);
	      
	      closed = 1;
	      sp_found = 1;
	      //TODO: ADD DISTANCE
	    }
	  }
	} 

	if (!closed) {
	  //SMITH-WATERMAN
	  if (gap_read <= 0 || gap_genome <= 0) {
	    seed_prev->genome_end   = seed_prev->genome_end - FLANK;
	    seed_next->genome_start = seed_next->genome_start + FLANK;
	    seed_prev->read_end   = seed_prev->read_end - FLANK;
	    seed_next->read_start = seed_next->read_start + FLANK;	    
	  }
	  size_t genome_start = seed_prev->genome_end + 1;
	  size_t genome_end   = seed_next->genome_start - 1;
	  int read_start      = seed_prev->read_end + 1;
	  int read_end        = seed_next->read_start - 1;

	  //assert(genome_end - genome_start + 1 < 2048);
	  if (genome_end - genome_start >= max_size) {
	    free(reference);
	    max_size = genome_end - genome_start + 1024;
	    reference = (char *)calloc(max_size, sizeof(char));
	  }
	  genome_read_sequence_by_chr_index(reference, 0, cal->chromosome_id - 1,
					    &genome_start, &genome_end, genome);
	  char query[2048];

	  assert(read_end - read_start + 1 < 2048);
	  assert(read_end - read_start + 1 > 0);

	  memcpy(query, &query_map[read_start],  read_end - read_start + 1);
	  query[read_end - read_start + 1] = '\0';

	  //printf("query : %s [%i-%i]\n", query, read_start, read_end);
	  //printf("ref   : %s [%lu-%lu]\n", reference, genome_start, genome_end);
	  seed_region_t *new_seed = seed_region_new(seed_prev->read_end + 1, seed_next->read_start - 1, 
						    seed_prev->genome_end + 1, seed_next->genome_start - 1, 
						    seed_prev->id + 1, 0, 0);
	  new_seed->fusion_left  = 1;
	  new_seed->fusion_right = 1;
	  linked_list_item_t *new_item = linked_list_item_new(new_seed);
	  list_item_prev->next = new_item;
	  new_item->prev = list_item_prev;
	  list_item->prev = new_item;
	  new_item->next = list_item;	  
	  cal->sr_list->size++;

	  sw_item = sw_item_new(SIMPLE_SW, i, 0, 0,
				cal, cal, NULL, 
				new_seed, new_seed,
				NULL);
	
	  //Insert item... and process if depth is full
	  sw_depth_insert(query, reference, sw_item,
			  sw_optarg, output,
			  avls_list, metaexons, 
			  sw_depth, genome, min_intron_size);	    

	  cal->num_targets++;
	}

	if (add_seed) {
	  seed_region_t *new_seed = seed_region_new(seed_prev->read_end + 1, seed_next->read_start - 1, 
						    seed_prev->genome_end + 1, seed_next->genome_start - 1, 
						    seed_prev->id + 1, 0, 0);
	  new_seed->fusion_left  = 1;
	  new_seed->fusion_right = 1;
	  linked_list_item_t *new_item = linked_list_item_new(new_seed);
	  list_item_prev->next = new_item;
	  new_item->prev = list_item_prev;
	  list_item->prev = new_item;
	  new_item->next = list_item;	  
	  cal->sr_list->size++;
	}

	if (!sp_found) {
	  seed_prev->fusion_right = 1;
	  seed_next->fusion_left  = 1;
	}

      }
      seed_prev      = seed_next;
      list_item_prev = list_item;
      list_item      = list_item->next;
    }

    cal->fill_gaps = 1;
    //printf("#### (%i)[%i:%lu-%lu] : %i ####\n", cal->strand, 
    //	   cal->chromosome_id, cal->start, cal->end, cal->num_targets);

    if (cal->num_targets == 0 && cal->fill_gaps) {
      //Close CAL
      merge_seeds_cal(cal);
    }
    //printf("END CAL %i\n", i); 
  }
  free(reference);
  //printf("EXIT LOOP\n");
}


meta_alignment_t *meta_alignment_cal_new(cal_t *cal) {
  meta_alignment_t *meta_alignment = meta_alignment_new();
  meta_alignment_insert_cal(cal, meta_alignment);

  return meta_alignment;
}

meta_alignment_t *meta_alignment_cals_new(array_list_t *cals_list) {
  meta_alignment_t *meta_alignment = meta_alignment_new();
  
  for (int i = 0; i < array_list_size(cals_list); i++) {
    cal_t *cal = array_list_get(i, cals_list);
    meta_alignment_insert_cal(cal, meta_alignment);
  }

  return meta_alignment;
}



//===============================================================================//


//============================= SMITH-WATERMAN SECTION ==========================//
info_sp_t *info_sp_new(size_t l_genome_start, size_t l_genome_end,
		       size_t r_genome_start, size_t r_genome_end) {

  info_sp_t *info_sp = (info_sp_t *)malloc(sizeof(info_sp_t));

  info_sp->l_genome_start = l_genome_start;
  info_sp->l_genome_end = l_genome_end;
  info_sp->r_genome_start = r_genome_start;
  info_sp->r_genome_end = r_genome_end;

  return info_sp;
}

void info_sp_free(info_sp_t *info_sp) {
  free(info_sp);
}

sw_item_t *sw_item_new(int type_sw, int read_id, 
		       int fusion_id, int cal_id,
		       cal_t *cal_prev, cal_t *cal_next, 
		       meta_alignment_t *meta_alignment, 
		       seed_region_t *seed_prev,
		       seed_region_t *seed_next,
		       void *info) {

  sw_item_t *sw_item = (sw_item_t *)malloc(sizeof(sw_item_t));
  
  sw_item->type_sw = type_sw;
  sw_item->read_id = read_id;
  sw_item->fusion_id = fusion_id;
  sw_item->cal_id = cal_id;
  sw_item->cal_prev = cal_prev;
  sw_item->cal_next = cal_next;
  sw_item->info = info;
  sw_item->meta_alignment = meta_alignment;
  sw_item->seed_prev = seed_prev;
  sw_item->seed_next = seed_next;

  return sw_item;

}

void sw_item_free(sw_item_t *sw_item) {
  free(sw_item);
}

void sw_depth_process(sw_optarg_t *sw_optarg, sw_multi_output_t *output, 
		      sw_depth_t *sw_depth, avls_list_t *avls_list,
		      metaexons_t *metaexons, int step, genome_t *genome,
		      int min_intron_size) {
  int distance, len;
  float norm_score;
  float match = sw_optarg->subst_matrix['A']['A'];
  if (sw_depth->depth == MAX_DEPTH || 
      (step == SW_FINAL && sw_depth->depth > 0)) {

    smith_waterman_mqmr(sw_depth->q, sw_depth->r, sw_depth->depth, sw_optarg, 1, output);

    //pthread_mutex_lock(&mutex_sp);
    //TOTAL_SW += sw_depth->depth;
    //pthread_mutex_unlock(&mutex_sp);

    for (int i = 0; i < sw_depth->depth; i++) {
      sw_item_t *sw_item = sw_depth->items[i];     
      cal_t *cal_prev = sw_item->cal_prev;
      cal_t *cal_next = sw_item->cal_next;
      
      //printf("-QUE: %s\n", sw_depth->q[i]);
      //printf("-REF: %s\n", sw_depth->r[i]);
      //printf("-QUE: %s(%i)\n", output->query_map_p[i], output->query_start_p[i]);
      //printf("-REF: %s(%i)\n", output->ref_map_p[i], output->ref_start_p[i]);
      
      if (sw_item->type_sw == EXTREM_SW_LEFT) {	
	norm_score = NORM_SCORE(output->score_p[i], strlen(sw_depth->q[i]), match);
	//printf("EXTREM SW LEFT SCORE %f\n", norm_score);
	cigar_code_t *cigar_code = NULL;
	if (norm_score >= 0.3) {
	  cigar_code = generate_cigar_code(output->query_map_p[i], output->ref_map_p[i],
					   strlen(output->ref_map_p[i]),
					   output->query_start_p[i], output->ref_start_p[i],
					   strlen(sw_depth->q[i]), strlen(sw_depth->r[i]),
					   &distance, FIRST_SW);
	  //printf("SW CIGAR %s\n", new_cigar_code_string(cigar_code));
	  //printf("....>%s\n", new_cigar_code_string(cigar_code));
	} 
	  //cigar_code_free(cigar_code);
	  //meta_alignment_insert_cigar(NULL, CIGAR_ANCHOR_RIGHT, sw_item->cal_id, sw_item->meta_alignment); 
	//}
	meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_RIGHT, sw_item->cal_id, sw_item->meta_alignment); 
      } else if (sw_item->type_sw == EXTREM_SW_RIGHT) {
	norm_score = NORM_SCORE(output->score_p[i], strlen(sw_depth->q[i]), match);
	//printf("EXTREM SW RIGHT SCORE %f\n", norm_score);
	cigar_code_t *cigar_code = NULL;
	if (norm_score >= 0.3) {
	  cigar_code = generate_cigar_code(output->query_map_p[i], output->ref_map_p[i],
					   strlen(output->ref_map_p[i]),
					   output->query_start_p[i], output->ref_start_p[i],
					   strlen(sw_depth->q[i]), strlen(sw_depth->r[i]),
					   &distance, LAST_SW);	
	  //printf("SW CIGAR %s\n", new_cigar_code_string(cigar_code));
	} 
	meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_LEFT, sw_item->cal_id, sw_item->meta_alignment); 
      } else if (sw_item->type_sw == SIMPLE_SW) {
	cigar_code_t *cigar_code = generate_cigar_code(output->query_map_p[i], output->ref_map_p[i],
						       strlen(output->ref_map_p[i]),
						       output->query_start_p[i], output->ref_start_p[i],
						       strlen(sw_depth->q[i]), strlen(sw_depth->r[i]),
						       &distance, MIDDLE_SW);
	cal_prev->num_targets--;
	seed_region_t *seed_prev = sw_item->seed_prev;
	seed_prev->info = cigar_code;
	
	//printf("SW CIGAR %s\n", new_cigar_code_string(cigar_code));
	//printf("(%i)[%i:%lu-%lu].....MERGE SEEDS %i??\n", cal_prev->strand, cal_prev->chromosome_id, 
	//     cal_prev->start, cal_prev->end, cal_prev->num_targets);
	if (cal_prev->num_targets == 0 && cal_prev->fill_gaps) {
	  //Close CAL
	  //printf("****MERGE SEEDS\n");
	  merge_seeds_cal(cal_prev);
	}
      } else {
	//fastq_read_t *read = array_list_get(sw_item->read_id, fq_batch);
	//if (read == NULL) { printf("@@@@@(%i)@ %s\n", sw_item->read_id, read->id); exit(-1); }
	info_sp_t *info_sp = sw_item->info;
	avl_node_t *node_avl_start, *node_avl_end;
	norm_score = NORM_SCORE(output->score_p[i], strlen(sw_depth->q[i]), match);
	//printf("QUE: %s(%i)\n", output->query_map_p[i], output->query_start_p[i]);
	//printf("REF: %s(%i)\n", output->ref_map_p[i], output->ref_start_p[i]);
	//printf("%f\n", norm_score);
	cigar_code_t *cigar_code;
	if (norm_score >= 0.6) {
	  cigar_code = generate_cigar_sw_output(output->query_map_p[i], 
						output->ref_map_p[i],
						info_sp->l_genome_start,
						info_sp->l_genome_end,
						info_sp->r_genome_start,
						info_sp->r_genome_end,
						sw_item->cal_prev->chromosome_id,
						sw_item->cal_prev->strand,
						output->query_start_p[i],
						output->ref_start_p[i],
						strlen(sw_depth->q[i]),
						strlen(sw_depth->r[i]),
						avls_list,
						&node_avl_start,
						&node_avl_end,
						genome, 
						min_intron_size);
	  //cigar_code = NULL;
	  //printf("CIGAR SW : %s\n", new_cigar_code_string(cigar_code));
	} else {
	  cigar_code = NULL;
	}
	info_sp_free(info_sp);
	if (cigar_code) {
	  int err_l = 0, err_r = 0;
	  if (node_avl_start && node_avl_end) {
	    
	    err_l = metaexon_insert(1, sw_item->cal_prev->chromosome_id - 1,
				    node_avl_start->position - 21, node_avl_start->position - 1, 40,
				    METAEXON_RIGHT_END, node_avl_start,
				    metaexons);
	    
	    err_r = metaexon_insert(1, sw_item->cal_prev->chromosome_id - 1,
				    node_avl_end->position + 1, node_avl_end->position + 21, 40,
				    METAEXON_LEFT_END, node_avl_end,
	    			    metaexons);
	    
	  }

	  if (err_l < 0 || err_r < 0) {
	    //printf("QUE: %s(%i)\n", output->query_map_p[i], output->query_start_p[i]);
	    //printf("REF: %s(%i)\n", output->ref_map_p[i], output->ref_start_p[i]);	    
	    exit(-1);
	  }
	}

	meta_alignment_insert_cigar(cigar_code, CIGAR_SW_MIDDLE, sw_item->cal_id, sw_item->meta_alignment);	  		
      }
      free(output->query_map_p[i]);
      free(output->ref_map_p[i]);
      output->query_map_p[i] = NULL;
      output->ref_map_p[i] = NULL;
      free(sw_depth->q[i]);
      free(sw_depth->r[i]);
      sw_item_free(sw_item);
    }
    sw_depth->depth = 0;
  }
}

void sw_depth_insert(char *query, char *reference, sw_item_t *sw_item, 
		     sw_optarg_t *sw_optarg, sw_multi_output_t *output, 
		     avls_list_t *avls_list, metaexons_t *metaexons, 
		     sw_depth_t *sw_depth, genome_t *genome, 
		     int min_intron_size) {

  sw_depth->q[sw_depth->depth]       = strdup(query);
  sw_depth->r[sw_depth->depth]       = strdup(reference);
  sw_depth->items[sw_depth->depth++] = sw_item;

  //printf("\tInsert SW depth %i\n", sw_depth->depth);
  //printf("QUE: %s\n", query);
  //printf("REF: %s\n", reference);

  sw_depth_process(sw_optarg, output, 
		   sw_depth, avls_list, metaexons, 
		   SW_NORMAL, genome, min_intron_size);

}
//===============================================================================//


info_sp_t* sw_reference_splice_junction(cal_t *cal_prev, cal_t *cal_next,
					char *query_map, genome_t *genome,
					char *q, char *r) {

  //printf("============= M O U N T    S M I T H - W A T E R M A N =============\n");  
  /*printf(" ==== CALS INFO ==== \n");
  cal_print(cal_prev);
  cal_print(cal_next);
  printf(" ==== CALS INFO END ==== \n");
  */
  seed_region_t *s_next, *s_prev;
  char reference[2048];
  char reference_prev[2048];
  char reference_next[2048];
  char query[2048];

  size_t genome_start, genome_end;
  int read_start, read_end, read_gap;
  size_t genome_start2, genome_end2;
  int seeds_nt;
  int flank = 30;
  int flank_left, flank_right;

  //Delete for debuging, detect splice junctions
  s_prev = linked_list_get_last(cal_prev->sr_list);    
  s_next = linked_list_get_first(cal_next->sr_list);

  read_start = s_prev->read_end;
  read_end = s_next->read_start;
  seeds_nt = read_end - read_start;

  //printf("SP coords %i - %i = %i\n", read_end, read_start, seeds_nt);

  if (seeds_nt >= 10) {
    //Extend to Right --> <-- Extend to Left
    genome_start = s_prev->genome_end;
    genome_end = s_prev->genome_end + seeds_nt - 1;
    genome_read_sequence_by_chr_index(reference_prev, 0, 
				      cal_prev->chromosome_id - 1, &genome_start, &genome_end, genome);

    genome_start2 = s_next->genome_start - seeds_nt;
    genome_end2 = s_next->genome_start - 1;
    genome_read_sequence_by_chr_index(reference_next, 0, 
				      cal_next->chromosome_id - 1, &genome_start2, &genome_end2, genome);

    memcpy(query, query_map + read_start, read_end - read_start);
    query[read_end - read_start]  = '\0';

    //CAll new function
    int dsp_e1, dsp_e2;
    int lim_err = 2;

    extend_by_mismatches(reference_prev, reference_next, query, 
			 0, strlen(reference_next) - 1, 
			 0, read_end - read_start - 1, lim_err,
			 &dsp_e1, &dsp_e2);
    //printf("dsp_e1=%i, dsp_e2=%i\n", dsp_e1, dsp_e2);

    cal_prev->end = cal_prev->end + dsp_e1;
    cal_next->start = cal_next->start - dsp_e2;

    s_prev->read_end = s_prev->read_end + dsp_e1;
    s_prev->genome_end = s_prev->genome_end + dsp_e1;

    s_next->read_start = s_next->read_start - dsp_e2;
    s_next->genome_start = s_next->genome_start - dsp_e2;

    if (s_prev->read_end > s_next->read_start) { 
      seeds_nt = s_prev->read_end - s_next->read_start;
      cal_prev->end -= seeds_nt;
      cal_next->start += seeds_nt;
      s_prev->read_end-= seeds_nt;
      s_next->read_start += seeds_nt;
    }
  }
    
  //printf("Read_start =: s_prev->read_end = %i - %i + 1 = %i\n", s_prev->read_end, flank, s_prev->read_end - flank + 1);
  
  int seed_size = s_prev->read_end - s_prev->read_start + 1;  
  flank_left = flank;
  if (flank > seed_size) {
    flank_left = seed_size;    
  } 
  read_start = s_prev->read_end - flank_left + 1;
  
  //printf("Read_end =: s_next->read_start = %i + %i - 1 = %i\n", s_next->read_start, flank, s_next->read_start + flank - 1);
  seed_size = s_next->read_end - s_next->read_start + 1;  
  flank_right = flank;
  if (flank > seed_size) {
    flank_right = seed_size;    
  } 
  read_end = s_next->read_start + flank_right - 1;
  

  //Extract and fusion Reference SW
  cal_prev->l_flank = flank_left;
  //cal_prev->r_flank = flank;
  genome_start = cal_prev->end - flank_left + 1;
  //printf("GENOME END %lu - %i + 1 = %lu\n", cal_prev->end, flank_left, genome_start);
  genome_end = cal_prev->end + flank - 1;
  genome_read_sequence_by_chr_index(reference_prev, 0, 
				    cal_prev->chromosome_id - 1, &genome_start, &genome_end, genome);
  //printf("1g[From %lu to %lu](%i): %s(%i)\n", genome_start, genome_end, genome_end - genome_start + 1, 
  //	 reference_prev, strlen(reference_prev));

  //cal_next->l_flank = flank;
  cal_next->r_flank = flank_right;
  genome_start2 = cal_next->start - flank;
  genome_end2 = cal_next->start + flank_right - 1;
  genome_read_sequence_by_chr_index(reference_next, 0, 
				    cal_next->chromosome_id - 1, &genome_start2, &genome_end2, genome);
  //printf("2g[From %lu to %lu](%i): %s(%i)\n", genome_start2, genome_end2, genome_end2 - genome_start2 + 1, 
  //	 reference_next, strlen(reference_next));

  strcat(reference_prev, reference_next);

  if (read_start > read_end) {
    int aux_start = read_end;
    read_start = read_end;
    read_end = aux_start;
  } else if (read_start == read_end) {
    LOG_FATAL("ERROR COORDS FUSION\n");
  }

  //if (read_start > read_end) { LOG_FATAL_F("READ COORDS ERROR %s\n", query_map); }
  //printf("Read %i-%i: %s\n", read_start, read_end, query_map);

  read_gap = read_end - read_start + 1;
  memcpy(query, &query_map[read_start], read_gap);
  query[read_gap] = '\0';

  //CALs Actualization flank
  //cal_prev->end -= flank;
  //cal_next->start += flank;
  //s_prev->read_end-= flank;
  //s_next->read_start += flank;
  //printf("\tread_start = %i, read_end = %i, genome_start = %lu, genome_end = %lu, flank_left=%i, flank_right=%i, query=%i =? read_gap=%i\n", 
  //	 read_start, read_end, genome_start, genome_end,  
  //	 flank_left, flank_right, strlen(query), read_gap);

  strcpy(q, query);
  strcpy(r, reference_prev);

  //printf("============= M O U N T    S M I T H - W A T E R M A N    E N D =============\n");

  return info_sp_new(genome_start, genome_end,
		     genome_start2, genome_end2);

}

//====================================================================================//

//avls_list_t *avls_list, 
//avl_node_t **node_avl_start,
//avl_node_t **node_avl_end,


//FOR SEARCH SEMI-CANNONICAL SPLICE JUNCTION WITHOUT SMITH-WATERMAN ALGORITHM
int search_simple_splice_junction_semi_cannonical(seed_region_t *s_prev, seed_region_t *s_next,
						  int chromosome_id, int strand, 
						  char *sequence, genome_t *genome, 
						  size_t *sp_start, size_t *sp_end,
						  int *sp_type,
						  int *distance) {

  assert(s_prev != NULL);
  assert(s_next != NULL);

  int seq_len = strlen(sequence);
  int read_start   = s_prev->read_end;
  int read_end     = s_next->read_start;
  int intron_size = 0;

  //printf("START READ START %i/%lu READ END %i/%lu\n", read_start, s_prev->genome_end, 
  //	 read_end, s_next->genome_start);

  size_t genome_start;
  size_t genome_end;
  
  const int FLANK = 10;
  const int SEQURITY_FLANK = 5;

  int gap_read = read_end - read_start - 1;
  if (gap_read == 0) {
    gap_read = -1;
  }
  
  int read_end_aux = s_prev->read_end;
  int read_start_aux = s_next->read_start;
  size_t genome_start_aux = s_next->genome_start;
  size_t genome_end_aux = s_prev->genome_end;

  LOG_DEBUG_F("SEARCH | %i-%i | %lu-%lu | \n", read_end_aux, read_start_aux, genome_end_aux, genome_start_aux);

  if (gap_read < 0) {
    //SEQURITY_FLANK = 5;
    gap_read = abs(gap_read) + SEQURITY_FLANK;
    read_end_aux     -= gap_read;
    read_start_aux   += gap_read;
    genome_end_aux   -= gap_read;
    genome_start_aux += gap_read;    
  } else {
    //SEQURITY_FLANK = gap_read + SEQURITY_FLANK;
    read_end_aux     -= SEQURITY_FLANK;
    read_start_aux   += SEQURITY_FLANK;
    genome_end_aux   -= SEQURITY_FLANK;
    genome_start_aux += SEQURITY_FLANK;
  }

  read_start = read_end_aux;
  read_end = read_start_aux;
  gap_read = read_end - read_start - 1;

  char left_exon[2048];
  char right_exon[2048];

  genome_start = genome_end_aux + 1;
  genome_end   = genome_end_aux + gap_read + FLANK;

  LOG_DEBUG_F("GAP READ %i - %i = %i\n", read_end, read_start, gap_read);
  LOG_DEBUG_F("SEQUENCE   : %s\n", sequence);

  genome_read_sequence_by_chr_index(left_exon, 0, 
				    chromosome_id - 1, 
				    &genome_start, &genome_end, genome);		  

  LOG_DEBUG_F("LEFT EXON  (%lu-%lu): %s\n", genome_start, genome_end, left_exon);

  genome_start = genome_start_aux - gap_read - FLANK;
  genome_end   = genome_start_aux - 1;
  
  genome_read_sequence_by_chr_index(right_exon, 0, 
				    chromosome_id - 1, 
				    &genome_start, &genome_end, genome);		  

  LOG_DEBUG_F("RIGHT EXON (%lu-%lu): %s\n", genome_start, genome_end, right_exon);

  size_t dsp_l, dsp_r, type;

  size_t breaks_starts[gap_read];
  int found_starts = 0, found_starts_semi = 0;

  size_t breaks_ends[gap_read];
  int found_ends = 0, found_ends_semi = 0;

  size_t type_starts[gap_read];
  size_t type_ends[gap_read];

  size_t breaks_starts_semi[gap_read];
  size_t breaks_ends_semi[gap_read];

  size_t type_starts_semi[gap_read];
  size_t type_ends_semi[gap_read];

  int c_s, c_e;

  int end_search;// = gap_read + SEQURITY_FLANK;

  if (strlen(left_exon) > strlen(right_exon)) { 
    end_search = strlen(right_exon);
  } else {
    end_search = strlen(left_exon);
  }
 
  //printf("search start!\n");
  // Search step by step (GT)/(AG) 
  for (c_s = 0, c_e = strlen(right_exon) - 1;
       c_s < end_search; c_s++, c_e--) {

    //Search Start Marks
    if (left_exon[c_s] == 'G' && left_exon[c_s + 1] == 'T') {
      type_starts[found_starts] = GT_AG_SPLICE;
      breaks_starts[found_starts++] = c_s;

      type_starts_semi[found_starts_semi] = GT_AT_SPLICE;
      breaks_starts_semi[found_starts_semi++] = c_s;
      LOG_DEBUG_F("S: FOUND GT (%i)\n", c_s);
    } else if (left_exon[c_s] == 'C' && left_exon[c_s + 1] == 'T') {
      type_starts[found_starts] = CT_AC_SPLICE;
      breaks_starts[found_starts++] = c_s;

      type_starts_semi[found_starts_semi] = CT_GC_SPLICE;
      breaks_starts_semi[found_starts_semi++] = c_s;
      LOG_DEBUG_F("S: FOUND CT (%i)\n", c_s);
    } else if (left_exon[c_s] == 'A' && left_exon[c_s + 1] == 'T') {
      LOG_DEBUG_F("S: FOUND AT (%i)\n", strlen(right_exon) - c_e - 1);
      type_starts_semi[found_starts_semi] = AT_AC_SPLICE;
      breaks_starts_semi[found_starts_semi++] = c_s;
    } else if (left_exon[c_s] == 'G' && left_exon[c_s + 1] == 'C') {
      LOG_DEBUG_F("S: FOUND GC (%i)\n", strlen(right_exon) - c_e - 1);
      type_starts_semi[found_starts_semi] = GC_AG_SPLICE;
      breaks_starts_semi[found_starts_semi++] = c_s;
    }

    //Search End Marks
    if (right_exon[c_e - 1] == 'A' && right_exon[c_e] == 'G') {
      type_ends[found_ends] = GT_AG_SPLICE;
      breaks_ends[found_ends++] = strlen(right_exon) - c_e - 1;

      type_ends_semi[found_ends_semi] = GC_AG_SPLICE;
      breaks_ends_semi[found_ends_semi++] = strlen(right_exon) - c_e - 1;      
      LOG_DEBUG_F("E: FOUND AG (%i)\n", strlen(right_exon) - c_e - 1);
    } else if (right_exon[c_e - 1] == 'A' && right_exon[c_e] == 'C') {
      type_ends[found_ends] = CT_AC_SPLICE;
      breaks_ends[found_ends++] = strlen(right_exon) - c_e - 1;

      type_ends_semi[found_ends_semi] = AT_AC_SPLICE;
      breaks_ends_semi[found_ends_semi++] = strlen(right_exon) - c_e - 1;      
      LOG_DEBUG_F("E: FOUND AC (%i)\n", strlen(right_exon) - c_e - 1);
    } else if (right_exon[c_e - 1] == 'A' && right_exon[c_e] == 'T') {
      type_ends_semi[found_ends_semi] = GT_AT_SPLICE;
      breaks_ends_semi[found_ends_semi++] = strlen(right_exon) - c_e - 1;
      LOG_DEBUG_F("E: FOUND AT (%i)\n", strlen(right_exon) - c_e - 1);
    } else if (right_exon[c_e - 1] == 'G' && right_exon[c_e] == 'C') {
      LOG_DEBUG_F("E: FOUND GC (%i)\n", strlen(right_exon) - c_e - 1);
      type_ends_semi[found_ends_semi] = CT_GC_SPLICE;
      breaks_ends_semi[found_ends_semi++] = strlen(right_exon) - c_e - 1;      
    }

  }

  //Not found any splice junction Cannonical
  if ((found_starts <= 0 || found_ends <= 0) && 
      (found_starts_semi <= 0 || found_ends_semi <= 0)) {
    //Not found any splice junction Semmi-Cannonical
    return 0;
  }   

  array_list_t *splice_junction = array_list_new(found_starts + found_ends, 
						 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  //array_list_t *splice_junction_1 = array_list_new(found_starts + found_ends, 
  //						   1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  //printf("FOUND STARTS %i, FOUND ENDS %i, FOUND STARTS SEMI %i, FOUND ENDS SEMI %i\n", found_starts, found_ends, found_starts_semi, 
  //	 found_ends_semi);


  //If found more than one break...Search cannonical
  for (int i = 0; i < found_starts; i++) {
    for (int j = 0; j < found_ends; j++) {
      //printf("%i(%i) + %i(%i)=%i\n", breaks_starts[i], type_starts[i],
      //     breaks_ends[j], type_ends[j], breaks_starts[i] + breaks_ends[j]);
      if (type_starts[i] == type_ends[j]) {
	int gap_break = breaks_starts[i] + breaks_ends[j];
	if (gap_break == gap_read) {
	  array_list_insert((void *)breaks_starts[i], splice_junction);
	  array_list_insert((void *)breaks_ends[j], splice_junction);
	  array_list_insert((void *)type_ends[j], splice_junction);
	}      
      }
    }
  }

  
  //If not found cannonical sp, search semmi-cannonical
  if (array_list_size(splice_junction) <= 0) {
    for (int i = 0; i < found_starts_semi; i++) {
      for (int j = 0; j < found_ends_semi; j++) {
	//printf("%i(%i) + %i(%i)=%i\n", breaks_starts_semi[i], type_starts_semi[i],
	//     breaks_ends_semi[j], type_ends_semi[j], breaks_starts_semi[i] + breaks_ends_semi[j]);
	if (type_starts_semi[i] == type_ends_semi[j]) {
	  int gap_break = breaks_starts_semi[i] + breaks_ends_semi[j];
	  if (gap_break == gap_read) {
	    array_list_insert((void *)breaks_starts_semi[i], splice_junction);
	    array_list_insert((void *)breaks_ends_semi[j], splice_junction);
	    array_list_insert((void *)type_ends_semi[j], splice_junction);
	  }      
	}
      }
    } 
  }


  //If not found splice-junctions exit.
  if (array_list_size(splice_junction) <= 0) {
    array_list_free(splice_junction, (void *)NULL);
    return 0;
  }


  //Calculating scores...
  int matches[array_list_size(splice_junction)];
  int mismatches[array_list_size(splice_junction)];
  float max_score = 0.0f;
  float score;
  int max_sp_pos;
  int sp_pos;
  int read_pos;
  int genome_pos;

  //printf("NUM SP (%i)\n", array_list_size(splice_junction)/3);
  for (int i = 0; i < array_list_size(splice_junction); i += 3) {
    int limit_left  = (size_t)array_list_get(i, splice_junction);
    int limit_right = (size_t)array_list_get(i + 1, splice_junction);      

    sp_pos = 0;
      
    matches[sp_pos] = 0;
    mismatches[sp_pos] = 0;

    read_pos = read_start + 1;
    for (int c_l = 0; c_l < limit_left; c_l++) {
      //printf("l: %c == %c\n", left_exon[c_l], sequence[c_l]);
      if (read_pos >= seq_len) { goto exit; }

      if (left_exon[c_l] == sequence[read_pos++]) {
	matches[sp_pos]++; 
      } else {
	mismatches[sp_pos]++; 
      }
    }

    read_pos = read_end - 1;
    genome_pos = strlen(right_exon) - 1;

    for (int c_r = 0; c_r < limit_right; c_r++) {
      //printf("r: %c == %c\n", right_exon[genome_pos], sequence[read_pos]);
      if (read_pos < 0) { goto exit; }

      if (right_exon[genome_pos--] == sequence[read_pos--]) {
	matches[sp_pos]++;
      } else {
	mismatches[sp_pos]++;
      }
    }

    score = matches[sp_pos]*0.5 - mismatches[sp_pos]*0.4;
    if (score > max_score) {
      max_sp_pos = sp_pos;
      max_score = score;
    }
   
    //printf("MATCHES %i / MISMATCHES %i\n", matches[sp_pos], mismatches[sp_pos]);
    sp_pos++;
  }

  dsp_l = (size_t)array_list_get((max_sp_pos * 3), splice_junction);
  dsp_r = (size_t)array_list_get((max_sp_pos * 3) + 1, splice_junction);
  type  = (size_t)array_list_get((max_sp_pos * 3) + 2, splice_junction);

  assert(type != 0);

  *distance = mismatches[max_sp_pos];

  int max_distance = 10;
  
  //printf("DISTANCE %i\n", *distance);

  if (type == CT_AC_SPLICE) {
    strand = 1;
  } else {
    strand = 0;
  }

  //assert(array_list_size(splice_junction) != 0);

  //TODO: CALCULATE DISTANCE 
  //printf("dsp_l = %i, dsp_r = %i, read_end_aux = %i, read_start_aux = %i\n", dsp_l, dsp_r, read_end_aux, read_start_aux);

  read_end_aux     += dsp_l;
  read_start_aux   -= dsp_r;
  genome_end_aux   += dsp_l;
  genome_start_aux -= dsp_r;
  

  if (*distance > max_distance || 
      (int)genome_start_aux - (int)genome_end_aux - 1 < 40) {
    intron_size = 0;
    goto exit;
  }

  s_prev->read_end     = read_end_aux;
  s_next->read_start   = read_start_aux;
  s_next->genome_start = genome_start_aux;
  s_prev->genome_end   = genome_end_aux;
  
  size_t start_splice = s_prev->genome_end + 1;
  size_t end_splice   = s_next->genome_start - 1;
  
  *sp_start = start_splice;
  *sp_end   = end_splice;
  *sp_type  = type;

  //printf("SP :=> [%i:%lu-%lu]\n", chromosome_id, start_splice, end_splice);
  intron_size = end_splice - start_splice + 1;
  
 exit:
  array_list_free(splice_junction, NULL);

  return intron_size;

}

//FOR SEARCH CANNONICAL SPLICE JUNCTION WITHOUT SMITH-WATERMAN ALGORITHM
int search_simple_splice_junction(seed_region_t *s_prev, seed_region_t *s_next,
				  int chromosome_id, int strand, 
				  char *sequence, genome_t *genome, 
				  size_t *sp_start, size_t *sp_end,
				  int *sp_type,
				  int *distance) {

  assert(s_prev != NULL);
  assert(s_next != NULL);

  int seq_len = strlen(sequence);
  int read_start   = s_prev->read_end;
  int read_end     = s_next->read_start;
  int intron_size = 0;

  //printf("START READ START %i/%lu READ END %i/%lu\n", read_start, s_prev->genome_end, 
  //	 read_end, s_next->genome_start);

  size_t genome_start;
  size_t genome_end;
  
  const int FLANK = 20;
  const int SEQURITY_FLANK = 5;

  int gap_read = read_end - read_start - 1;
  if (gap_read == 0) {
    gap_read = -1;
  }
  
  int read_end_aux = s_prev->read_end;
  int read_start_aux = s_next->read_start;
  size_t genome_start_aux = s_next->genome_start;
  size_t genome_end_aux = s_prev->genome_end;

  LOG_DEBUG_F("SEARCH | %i-%i | %lu-%lu | \n", read_end_aux, read_start_aux, genome_end_aux, genome_start_aux);
  if (gap_read < 0) {
    gap_read = abs(gap_read) + 5;
    read_end_aux     -= gap_read;
    read_start_aux   += gap_read;
    genome_end_aux   -= gap_read;
    genome_start_aux += gap_read;    
  } else {
    read_end_aux     -= SEQURITY_FLANK;
    read_start_aux   += SEQURITY_FLANK;
    genome_end_aux   -= SEQURITY_FLANK;
    genome_start_aux += SEQURITY_FLANK;
  }

  
  read_start = read_end_aux;
  read_end = read_start_aux;
  
  gap_read = read_end - read_start - 1;
  //printf("%i - %i - 1 = %i\n", read_end, read_start, gap_read);

  char left_exon[2048];
  char right_exon[2048];

  genome_start = genome_end_aux + 1;
  genome_end   = genome_end_aux + gap_read + FLANK;

  LOG_DEBUG_F("GAP READ %i - %i = %i\n", read_end, read_start, gap_read);
  LOG_DEBUG_F("SEQUENCE   : %s\n", sequence);

  LOG_DEBUG_F("GAP READ %i - %i = %i\n", read_end, read_start, gap_read);
  genome_read_sequence_by_chr_index(left_exon, 0, 
				    chromosome_id - 1, 
				    &genome_start, &genome_end, genome);		  

  LOG_DEBUG_F("LEFT EXON  (%lu-%lu): %s\n", genome_start, genome_end, left_exon);

  genome_start = genome_start_aux - gap_read - FLANK;
  genome_end   = genome_start_aux - 1;
  
  genome_read_sequence_by_chr_index(right_exon, 0, 
				    chromosome_id - 1, 
				    &genome_start, &genome_end, genome);		  

  LOG_DEBUG_F("RIGHT EXON (%lu-%lu): %s\n", genome_start, genome_end, right_exon);

  size_t dsp_l, dsp_r, type;

  size_t breaks_starts[gap_read];
  int found_starts = 0;

  size_t breaks_ends[gap_read];
  int found_ends = 0;

  size_t type_starts[gap_read];
  size_t type_ends[gap_read];
  int c_s, c_e;
  int end_search = gap_read + SEQURITY_FLANK;

  //printf("search start!\n");
  // Search step by step (GT)/(AG) 
  for (c_s = 0, c_e = strlen(right_exon) - 1;
       c_s < end_search; c_s++, c_e--) {
    if (left_exon[c_s] == 'G' && left_exon[c_s + 1] == 'T') {
      type_starts[found_starts] = GT_AG_SPLICE;
      breaks_starts[found_starts++] = c_s;
      LOG_DEBUG_F("S.FOUND GT (%i)\n", c_s);
    } else if (left_exon[c_s] == 'C' && left_exon[c_s + 1] == 'T') {
      type_starts[found_starts] = CT_AC_SPLICE;
      breaks_starts[found_starts++] = c_s;
      LOG_DEBUG_F("S.FOUND CT (%i)\n", c_s);
    }

    if (right_exon[c_e] == 'G' && right_exon[c_e - 1] == 'A') {
      type_ends[found_ends] = GT_AG_SPLICE;
      breaks_ends[found_ends++] = strlen(right_exon) - c_e - 1;
      LOG_DEBUG_F("E.FOUND AG (%i)\n", strlen(right_exon) - c_e - 1);
    } else if (right_exon[c_e] == 'C' && right_exon[c_e - 1] == 'A') {
      type_ends[found_ends] = CT_AC_SPLICE;
      breaks_ends[found_ends++] = strlen(right_exon) - c_e - 1;
      LOG_DEBUG_F("E.FOUND AC (%i)\n", strlen(right_exon) - c_e - 1);
    }
  }

  //Not found any splice junction
  if (found_starts == 0 || found_ends == 0) {
    return 0;
  } 

  array_list_t *splice_junction = array_list_new(found_starts + found_ends, 
						 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  //array_list_t *splice_junction_1 = array_list_new(found_starts + found_ends, 
  //						   1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  //printf("FOUND STARTS %i, FOUND ENDS %i\n", found_starts, found_ends);
  //printf("Left  exon: %s\n", left_exon);
  //printf("Right exon: %s\n", right_exon);

  //If found more than one break...
  for (int i = 0; i < found_starts; i++) {
    for (int j = 0; j < found_ends; j++) {
      //printf("%i(%i) + %i(%i)=%i\n", breaks_starts[i], type_starts[i],
      //breaks_ends[j], type_ends[j], breaks_starts[i] + breaks_ends[j]);
      if (type_starts[i] == type_ends[j]) {
	int gap_break = breaks_starts[i] + breaks_ends[j];
	//printf("::(%i, %i), ¿%i - %i?\n", breaks_starts[i], breaks_ends[j], gap_break, gap_read);
	if (gap_break == gap_read) {	  
	  //printf("END:%i-START:%i\n", type_ends[j], type_starts[i]);
	  array_list_insert((void *)breaks_starts[i], splice_junction);
	  array_list_insert((void *)breaks_ends[j], splice_junction);
	  array_list_insert((void *)type_ends[j], splice_junction);
	}      
      }
    }
  }

  if (array_list_size(splice_junction) <= 0) {
    array_list_free(splice_junction, (void *)NULL);
    return 0;
  }

  //Calculating scores...
  int matches[array_list_size(splice_junction)];
  int mismatches[array_list_size(splice_junction)];
  float max_score = 0.0f;
  float score = 0.0f;
  int max_sp_pos = 0;
  int sp_pos;
  int read_pos;
  int genome_pos;

  sp_pos = 0;
  //printf("NUM SP (%i)\n", array_list_size(splice_junction)/3);
  for (int i = 0; i < array_list_size(splice_junction); i += 3) {
    int limit_left  = (size_t)array_list_get(i, splice_junction);
    int limit_right = (size_t)array_list_get(i + 1, splice_junction);      

    //printf("%i vs %i\n", limit_left, limit_right);
      
    matches[sp_pos] = 0;
    mismatches[sp_pos] = 0;

    read_pos = read_start + 1;
    for (int c_l = 0; c_l < limit_left; c_l++) {
      //printf("l: %c == %c\n", left_exon[c_l], sequence[c_l]);
      if (read_pos >= seq_len) { goto exit; }

      if (left_exon[c_l] == sequence[read_pos++]) {
	matches[sp_pos]++; 
      } else {
	mismatches[sp_pos]++; 
      }

    }

    read_pos = read_end - 1;
    genome_pos = strlen(right_exon) - 1;

    for (int c_r = 0; c_r < limit_right; c_r++) {
      //printf("r: %c == %c\n", right_exon[genome_pos], sequence[read_pos]);
      if (read_pos < 0) { goto exit; }

      if (right_exon[genome_pos--] == sequence[read_pos--]) {
	matches[sp_pos]++;
      } else {
	mismatches[sp_pos]++;
      }

    }

    //printf("matches=%i - mismatches=%i\n", matches[sp_pos], mismatches[sp_pos]);

    score = (float)matches[sp_pos]*0.5 - (float)mismatches[sp_pos]*0.4;

    //printf("%f - \n", score);

    if (score > max_score) {
      //printf("SP-YYPE: %i\n", type_ends[sp_pos]);
      max_sp_pos = sp_pos;
      max_score = score;
    }
   
    //printf("MATCHES %i / MISMATCHES %i\n", matches[sp_pos], mismatches[sp_pos]);
    sp_pos++;

  }

  if (max_score < 0.0) { 
    //printf("Exit\n");
    intron_size = 0;
    goto exit;
  }

  dsp_l = (size_t)array_list_get((max_sp_pos * 3), splice_junction);
  dsp_r = (size_t)array_list_get((max_sp_pos * 3) + 1, splice_junction);
  type  = (size_t)array_list_get((max_sp_pos * 3) + 2, splice_junction);

  *distance = mismatches[max_sp_pos];

  int max_distance = 10;
  
  //printf("DISTANCE %i\n", *distance);

  if (type == CT_AC_SPLICE) {
    strand = 1;
  } else {
    strand = 0;
  }

  //assert(array_list_size(splice_junction) != 0);

  //TODO: CALCULATE DISTANCE 
  //printf("dsp_l = %i, dsp_r = %i, read_end_aux = %i, read_start_aux = %i\n", dsp_l, dsp_r, read_end_aux, read_start_aux);

  read_end_aux     += dsp_l;
  read_start_aux   -= dsp_r;
  genome_end_aux   += dsp_l;
  genome_start_aux -= dsp_r;
  

  if (*distance > max_distance || 
      (int)genome_start_aux - (int)genome_end_aux - 1 < 40) {
    intron_size = 0;
    goto exit;
  }

  s_prev->read_end     = read_end_aux;
  s_next->read_start   = read_start_aux;
  s_next->genome_start = genome_start_aux;
  s_prev->genome_end   = genome_end_aux;
  
  size_t start_splice = s_prev->genome_end + 1;
  size_t end_splice   = s_next->genome_start - 1;
  
  *sp_start = start_splice;
  *sp_end   = end_splice;
  *sp_type  = type;

  //printf("SP :=> [%i:%lu-%lu]\n", chromosome_id, start_splice, end_splice);
  intron_size = end_splice - start_splice + 1;
  
 exit:
  array_list_free(splice_junction, NULL);

  return intron_size;

}

cigar_code_t *search_left_single_anchor(int gap_close, 
					cal_t *cal,
					int filter_pos, 
					array_list_t *right_breaks,
					char *query_map,
					metaexons_t *metaexons,
					genome_t *genome,
					avls_list_t *avls_list) {
  //return NULL;  
  //cal_print(cal);
  cigar_code_t *cigar_code = NULL;
  int max_nt;
  avl_node_t *node_start_prev, *node_start_next;
  int dist, final_dist = 0;
  int abort = 0, map = 0;
 
  array_list_t *final_positions = array_list_new(50, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *starts_targets  = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *final_starts    = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
 
  metaexon_t *final_metaexon;
  size_t final_pos;
  char reference[2048];
  size_t genome_start, genome_end;
  avl_node_t *node_start;
  int max_dist;
  int s;
  seed_region_t *s_prev = linked_list_get_last(cal->sr_list);

  int read_pos = s_prev->read_end + 1;

  size_t first_cal_end = cal->end;
  int cal_strand = cal->strand;
  int cal_chromosome_id = cal->chromosome_id;

  genome_start = first_cal_end + 1;
  //printf("IN PARAMETERS: read_pos->%i, first_cal_end->%lu, gap_close->%i, genome_start->%lu\n", read_pos, first_cal_end, gap_close, genome_start );
  //array_list_insert(genome_start, final_positions);
  //==== 1st-Select the correct start ====//
  for (int s_0 = 0; s_0 < array_list_size(right_breaks); s_0++) {
    node_start = array_list_get(s_0, right_breaks);
    if (first_cal_end - 10 <= node_start->position) {
      //LOG_DEBUG_F("\tNode start %lu\n", node_start->position);
      array_list_insert(node_start, starts_targets);
    }
  }
  
  while (gap_close > 0) {
    //1. Order starts
    int num_targets = array_list_size(starts_targets); 
    if (!num_targets) { 
      //printf("Not Targets\n");
      //Exonic read?      
      genome_end = genome_start + gap_close + 1;

      //printf("%i + %i (%i) < %i\n", read_pos, gap_close,
      //	     read_pos + gap_close, strlen(query_map));
      //printf("%lu - %lu\n", genome_start, genome_end);

      assert((genome_end - genome_start) + 5 < 2048);
      
      //printf("::: %i+%i=%i <= %i\n", read_pos, gap_close, read_pos + gap_close, strlen(query_map));
      assert((read_pos + gap_close) <= strlen(query_map));

      genome_read_sequence_by_chr_index(reference, 0, 
					cal_chromosome_id - 1, 
					&genome_start, &genome_end, genome);
      int reference_len = strlen(reference) - 1;
      dist = 0;
      //printf("Not Splice Exonic Read...cheking\n");
      int c, lim_ref = gap_close;
      for (c  = 0; c < gap_close; c++) { 
	//printf("\t\t[%c vs %c]\n", query_map[c], reference[c]);		
	if (query_map[c + read_pos] != reference[c]) { 
	  dist++;
	}
      }

      //max_dist = lim_ref <= 5 ? 2 : lim_ref / 4 + 1;
      max_dist = lim_ref <= 5 ? 2 : lim_ref / 4;
      //printf("%i < %i\n", dist, max_dist);
      if (dist < max_dist) {
	final_dist += dist;
	map = 1;
	array_list_insert((void *)genome_start + gap_close - 1, final_positions);
	//printf(" Insert 0) %lu\n", genome_start + gap_close - 1);
      }
      break; 
    }
    for (int s0 = 0; s0 < num_targets - 1; s0++) {
      node_start_prev = array_list_get(s0, starts_targets);
      for (int s1 = s0 + 1; s1 <  num_targets; s1++) {
	node_start_next = array_list_get(s1, starts_targets);
	if (node_start_next->position < node_start_prev->position) {
	  array_list_swap(s0, s1, starts_targets);
	}
      } 
    }

    //2. Map the gap to the genome and Select the correct start node
    node_start = array_list_get(num_targets - 1, starts_targets);
    genome_end = node_start->position;
    int lim_ref;
    lim_ref = genome_end - genome_start;
    //printf("genome_end = %lu, genome_start = %lu, lim_ref = %i\n", genome_end, genome_start, lim_ref);

    if (lim_ref > gap_close) {
      lim_ref = gap_close;
      genome_end = genome_start + lim_ref + 1;
    }

    if (lim_ref < 0) {
      //For first step, genome_start = first_cal->end
      int dsp = abs(lim_ref);
      gap_close += dsp;
      read_pos  -= dsp;
      array_list_insert((void *)genome_end - 1, final_positions);
      //printf("Read pos %i\n ", read_pos);
      //printf("1):::::::::::: INSERT FINAL POSITIONS: %i, (%i)\n", genome_end, array_list_size(final_positions));
      //printf(":::--::::GENOME GAP %i\n", genome_end);
      //printf("RECALCULATING GAP AND READ POS (gap_close, read_pos)(%i, %i)\n", gap_close, read_pos);
    } else {
      assert((genome_end - genome_start) + 5 < 2048);

      genome_read_sequence_by_chr_index(reference, 0, 
					cal_chromosome_id - 1, 
					&genome_start, &genome_end, genome);
      //printf("(%i)after extraction genome_start = %lu, genome_end = %lu\n", cal_chromosome_id, genome_start, genome_end);

      //printf("(Read Pos %i) Reference START_SP [[[START_SP]]]----[END_SP] Targets %i: %s\n", read_pos, num_targets, reference);

      int t, c;
      for (t = 0; t < num_targets; t++) {
	node_start = array_list_get(t, starts_targets);
	//printf("Calculating lim_ref = %lu - %lu\n", node_start->position, genome_start);
	lim_ref = node_start->position - genome_start;
	if (lim_ref > gap_close) { lim_ref = gap_close; }
	dist = 0;

	//assert((read_pos + lim_ref) < strlen(query_map));
	if ((read_pos + lim_ref) > strlen(query_map)) {
	  LOG_FATAL_F("%i + %i (%i) < %i\n", read_pos, lim_ref,
		      read_pos + lim_ref, strlen(query_map));
	}

	//printf("\tTravel to %lu, gap_close=%i, lim_ref=%i\n", node_start->position, gap_close, lim_ref);
	for (c  = 0; c < lim_ref; c++) { 
	  //printf("\t\t[%c vs %c]\n", query_map[read_pos + c], reference[c]);
	  if (read_pos + c >= strlen(query_map)) { printf("Overfflow 1");exit(-1); }
	  if (c >= strlen(reference)) { printf("Overfflow 1"); exit(-1); }

	  if (query_map[read_pos + c] != reference[c]) { 
	    dist++;
	  }
	}
	max_dist = lim_ref <= 5 ? 2 : lim_ref / 4;
	//printf("%i < %i\n", dist, max_dist);
	if (dist < max_dist) { 
	  //printf("\t\tMAX DISTANCE GOOD %i!\n", dist);
	  //max_dist += dist;
	  final_dist += dist;
	  if (lim_ref == gap_close) { t++; break; }
	} else { 
	  //printf("\t\tExit LOOP %i\n", t - 1);
	  break;
	}
      }
      t--;
      if (t < 0) { break; }
      else {
	node_start = array_list_get(t, starts_targets);
	lim_ref = node_start->position - genome_start;
	if (lim_ref > gap_close) { lim_ref = gap_close; }
	gap_close -= lim_ref;
	read_pos += lim_ref;	
	//printf("GAP CLOSE (%i) READ POS (%i)\n", gap_close, read_pos);
	if (gap_close <= 0) {
	  //Final Map
	  map = 1;
	  array_list_insert((void *)genome_start + c - 1, final_positions);
	  //printf("3):::::::::::: INSERT FINAL POSITIONS: %i (%i)\n", genome_start + c - 1, array_list_size(final_positions));
	  break;
	} else {
	  array_list_insert((void *)node_start->position - 1, final_positions);
	  //printf("2):::::::::::: INSERT FINAL POSITIONS: %i (%i)\n", node_start->position - 1, array_list_size(final_positions));
	}
      }
    }
	    
    //3.Select the correct start
    int pos;
    array_list_t *ends_list = ((start_data_t *)node_start->data)->list_ends;

    array_list_clear(final_starts, (void *)NULL);
    //Filter ends before cal_next->start
    for (int e = 0; e < array_list_size(ends_list); e++) {
      splice_end_t *splice_end = array_list_get(e, ends_list);
      size_t start = splice_end->end;
      if (filter_pos) {
	if (start < filter_pos) {
	  array_list_insert((void *)start, final_starts);
	}
      } else {
	array_list_insert((void *)start, final_starts);
      }
    }
    if (array_list_size(final_starts) >= 1) {
      //Select best start
      char s_reference[array_list_size(final_starts)][2048];
      //Making reference...
      for (int s = 0; s < array_list_size(final_starts); s++) {
	genome_start = (size_t)array_list_get(s, final_starts) + 1;
	genome_end = genome_start + gap_close + 1;

	assert((genome_end - genome_start) + 5 < 2048);
	//printf("%lu - %lu\n", genome_start, genome_end);
	genome_read_sequence_by_chr_index(s_reference[s], 0, 
					  cal_chromosome_id - 1, 
					  &genome_start, &genome_end, genome); 	      
	//printf("Reference END_SP [START_SP]----[[[END_SP]]]  [%i:%lu-%lu](%i): %s\n", cal_chromosome_id, genome_start, genome_end, 
	//s, s_reference[s]);
      } 

      //printf("GAP CLOSE (%i) READ POS (%i)\n", gap_close, read_pos);
      //Select the best start
      int anchor_nt = 20; 
      if (anchor_nt > gap_close) { 
	anchor_nt = gap_close;
      }

      int *ref_matches = (int *)calloc(array_list_size(final_starts), sizeof(int));
      int *ref_mismatches = (int *)calloc(array_list_size(final_starts), sizeof(int));
      //printf("anchor nt = %i, array_list_size(final_starts)=%i\n", anchor_nt, array_list_size(final_starts));
      for (int c = 0; c < anchor_nt; c++) {
	for (int s = 0; s < array_list_size(final_starts); s++) {
	  //assert(read_pos + c <= strlen(query_map));
	  if (read_pos + c > strlen(query_map)) { 
	    printf("ERROR: %i+%i=%i  vs %i", read_pos, c, read_pos + c, strlen(query_map)); cal_print(cal); 
	    goto free;
	    //exit(-1); 
	  }

	  assert(c <= strlen(s_reference[s]));
	  //printf("\t(%i )[%c vs %c]\n", s, query_map[read_pos + c], s_reference[s][c]);
	  if (query_map[read_pos + c] != s_reference[s][c]) {
	    ref_mismatches[s]++;
	  } else {
	    ref_matches[s]++;
	  }
	}
      }
      pos = 0;
      float score, max_score = 0.0;
      for (int s = 0; s < array_list_size(final_starts); s++) {
	score = ref_matches[s]*0.5 - ref_mismatches[s]*0.4;
	if (score > max_score) { 
	  max_score = score;
	  pos = s;
	}	
      }		
      max_dist = anchor_nt <= 5 ? 2 : anchor_nt / 4 + 1;
      //printf("%i > %i\n", ref_mismatches[pos], max_dist);
      if (ref_mismatches[pos] > max_dist) {
	//printf("ohhhhhhhhhhhhh!\n");
	free(ref_matches);
	free(ref_mismatches);
	break;
      }

      final_dist += ref_mismatches[pos]; 
      final_pos = (size_t)array_list_get(pos, final_starts);
      array_list_insert((void *)final_pos, final_positions);
      //printf("4):::::::::::: INSERT FINAL POSITIONS: %i (%i)\n", final_pos, array_list_size(final_positions));
      genome_start = final_pos + anchor_nt + 1;
      gap_close -= anchor_nt;
      read_pos += anchor_nt;
      //printf("FINAL MATCHES (%i) VS FINAL MISMATCHES (%i)\n", ref_matches[s], ref_mismatches[s]);
      //printf("GAP CLOSE (%i) READ POS (%i)\n", gap_close, read_pos);

      free(ref_mismatches);
      free(ref_matches);

      if (anchor_nt < 20) {
	array_list_insert((void *)final_pos + anchor_nt, final_positions);
	//printf("5):::::::::::: INSERT FINAL POSITIONS: %i (%i)\n", final_pos + anchor_nt, array_list_size(final_positions));
	map = 1;
	break; 
      } else {
	if (metaexon_search(cal_strand, cal_chromosome_id - 1, 
			    final_pos, final_pos + anchor_nt, 
			    &final_metaexon, metaexons)) {
	  if (final_metaexon) {
	    if (final_metaexon->right_closed) {
	      array_list_clear(starts_targets, (void *)NULL);
	      //printf("FOUND!!\n");
	      for (int s_0 = 0; s_0 < array_list_size(final_metaexon->right_breaks); s_0++) {
		node_start = array_list_get(s_0, final_metaexon->right_breaks);
		if (final_pos <= node_start->position) {
		  array_list_insert(node_start, starts_targets);		    
		}
	      }
	    } else {
	      //We can close the gap
	      int exon_length = final_metaexon->end - final_metaexon->start;
	      dist = 0;
	      if (exon_length >= gap_close) {
		genome_end = genome_start + gap_close + 1;
		assert((genome_end - genome_start) + 5 < 2048);
		assert((read_pos + gap_close) <= strlen(query_map));
		genome_read_sequence_by_chr_index(reference, 0, 
						  cal_chromosome_id - 1, 
						  &genome_start, &genome_end, genome);		
		//printf("CLOSE GAP: %s\n", reference);
		for (int c  = 0; c < gap_close; c++) { 
		  //printf("\t[%c vs %c]\n", query_map[read_pos + c], reference[c]);
		  if (query_map[read_pos + c] != reference[c]) { 
		    dist++;
		  }
		}
		if (dist < 5) {
		  array_list_insert((void *)final_pos + gap_close + anchor_nt, final_positions);
		  //printf("6):::::::::::: INSERT FINAL POSITIONS: %i (%i)\n", final_pos + gap_close + anchor_nt,
		  //	 array_list_size(final_positions));
		  final_dist += dist; 
		  map = 1;
		}
	      }
	      break;
	    }
	  }
	} else {
	  //printf("NOT  FOUND\n");
	  break;
	}
      }
    } else {
      break;
    }
  }


  if (map) {
    //printf("read map!!\n");
    size_t pos_prev = cal->start, pos_next;
    cigar_op_t *op;
    cigar_code = cigar_code_new();
    cigar_code->distance = final_dist;

    pos_prev = cal->end + 1;
    //printf("FINAL POSITION %i\n", array_list_size(final_positions));
    for (int sp = 0; sp < array_list_size(final_positions); sp++) {
      pos_next = (size_t)array_list_get(sp, final_positions);
      if (sp % 2 == 0) {		
	//printf("%lu - %lu(%i/%i) = %i\n", pos_prev, pos_next, sp, array_list_size(final_positions), pos_next - pos_prev + 1);
	op = cigar_op_new(pos_next - pos_prev + 1, 'M');
      } else {
	size_t start_sp = pos_prev;
	size_t end_sp = pos_next;
	int aux = end_sp - start_sp + 1;

	op = cigar_op_new(aux, 'N');

	char nt_end[5], nt_start[5];
	int type;

	//Report splice junction
	size_t g_start = start_sp;
	size_t g_end   = start_sp + 1;
	genome_read_sequence_by_chr_index(nt_start, cal->strand, cal->chromosome_id - 1, &g_start, &g_end, genome);

	g_end   = end_sp;
	g_start = end_sp - 1;
	genome_read_sequence_by_chr_index(nt_end, cal->strand, cal->chromosome_id - 1, &g_start, &g_end, genome);
		
	type = splice_junction_type(nt_start[0], nt_start[1], nt_end[0], nt_end[1]);
	int splice_strand;
	avl_node_t *avl_node_start, *avl_node_end;
	
	splice_strand = 0;
	if (type == CT_AC_SPLICE || type == GT_AT_SPLICE || type == CT_GC_SPLICE ) {
	  splice_strand = 1;
	}

	//printf("L:(%i)[%i:%lu-%lu] : %c%c vs %c%c\n", cal->strand, cal->chromosome_id, start_sp, end_sp, nt_start[0], nt_start[1], nt_end[0], nt_end[1]);

	//####-1
	char str_sp_type[10];
	if (type == NOT_SPLICE) {
	  type = UNKNOWN_SPLICE;
	  sprintf(str_sp_type, "%c%c-%c%c\0", nt_start[0], nt_start[1], nt_end[0], nt_end[1]);
	}
	
	allocate_start_node(cal->chromosome_id - 1,
			    splice_strand,
			    start_sp,
			    end_sp,
			    start_sp,
			    end_sp,
			    FROM_READ,
			    type,
			    str_sp_type, 
			    &avl_node_start,
			    &avl_node_end, 
			    avls_list);	
      }

      //printf("\tADD %i%c\n", op->number, op->name);
      cigar_code_append_op(op, cigar_code);
      pos_prev = pos_next + 1;
    }
  }

  //printf("Final Cigar %s\n", new_cigar_code_string(cigar_code));
 free:
  array_list_free(final_positions, NULL);
  array_list_free(starts_targets, NULL);
  array_list_free(final_starts, NULL);

  //if (cigar_code == NULL) { exit(-1); }

  return cigar_code;

}

cigar_code_t *search_right_single_anchor(int gap_close, 
					 cal_t *cal,
					 int filter_pos, 
					 array_list_t *left_breaks,
					 char *query_map, metaexons_t *metaexons,
					 genome_t *genome,
					 avls_list_t *avls_list) {
  //return NULL;
  
  int max_nt;
  avl_node_t *node_end_prev, *node_end_next;
  int dist, final_dist = 0;
  int abort = 0, map = 0;
 
  array_list_t *final_positions = array_list_new(50, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *ends_targets = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *final_ends   = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
 
  metaexon_t *final_metaexon;
  size_t final_pos;
  char reference[2048];
  size_t genome_start, genome_end;
  avl_node_t *node_end;
  int max_dist;
  int s;
  size_t reference_len;

  seed_region_t *s_prev = linked_list_get_first(cal->sr_list);
  int read_pos = s_prev->read_start - 1; 
  size_t last_cal_start =  cal->start;
  int cal_strand = cal->strand;
  int cal_chromosome_id = cal->chromosome_id;

  cigar_code_t *cigar_code = NULL;
  
  //printf(" =============:> IN PARAMETERS: read_pos->%i, first_cal_end->%lu, gap_close->%i\n", read_pos, last_cal_start - 1, gap_close );
  /*
  printf("=========================================()====================================\n");
  linked_list_t *list_meta = metaexons->metaexons_list[18];
  linked_list_item_t *item_aux = list_meta->first;
  while (item_aux != NULL) {
    metaexon_t *meta_aux= item_aux->item;
    if (meta_aux->start >= 17510000 && meta_aux->start <= 17528000) {
      printf(" - %i-%i - ", meta_aux->start, meta_aux->end);
    }

    item_aux = item_aux->next;
  }  
  printf("\n=========================================()====================================\n");
  */

  genome_end = last_cal_start -  1;
  //array_list_insert(genome_start, final_positions);
  //==== 1st-Select the correct start ====//
  for (int s_0 = 0; s_0 < array_list_size(left_breaks); s_0++) {
    node_end = array_list_get(s_0, left_breaks);
    //printf("<XXX [%i >= %i] XXX>\n", last_cal_start + 10, node_end->position);
    if (last_cal_start + 10 >= node_end->position) {
      //printf("\tNode end %lu\n", node_end->position);
      array_list_insert(node_end, ends_targets);
    }
  }
  
  while (gap_close > 0) {
    //printf("Start loop...\n");
    //1. Order starts
    int num_targets = array_list_size(ends_targets); 
    if (!num_targets) { 
      //Exonic read?
      genome_start = genome_end - gap_close - 1;
      assert((genome_end - genome_start) + 2 < 2048);
      //assert((read_pos + gap_close) <= strlen(query_map));

      genome_read_sequence_by_chr_index(reference, 0, 
					cal_chromosome_id - 1, 
					&genome_start, &genome_end, genome);
      reference_len = strlen(reference) - 1;
      dist = 0;
      //printf("Not Splice Exonic Read...cheking\n");
     
      int c, lim_ref = gap_close;
      for (c  = 0; c < gap_close; c++) { 
	//printf("\t\t[%c vs %c]\n", query_map[read_pos - c], reference[reference_len - c]);		
	assert(read_pos - c >= 0);
	if (query_map[read_pos - c] != reference[reference_len - c]) { 
	  dist++;
	}
      }
      //max_dist = lim_ref <= 5 ? 2 : lim_ref / 4 + 1;
      max_dist = lim_ref <= 5 ? 2 : lim_ref / 4;
      //printf("%i < %i\n", dist, max_dist);
      if (dist <= max_dist) {
	final_dist += dist;
	map = 1;
	array_list_insert((void *)genome_end - gap_close + 1, final_positions);
	//printf(" Insert 0) %lu\n", genome_end - gap_close + 1);
      }
      //printf("Not targets\n");       
      break;
    }

    for (int s0 = 0; s0 < num_targets - 1; s0++) {
      node_end_prev = array_list_get(s0, ends_targets);
      for (int s1 = s0 + 1; s1 <  num_targets; s1++) {
	node_end_next = array_list_get(s1, ends_targets);
	if (node_end_next->position > node_end_prev->position) {
	  array_list_swap(s0, s1, ends_targets);
	}
      } 
    }

    //2. Map the gap to the genome and Select the correct start node
    node_end = array_list_get(num_targets - 1, ends_targets);
    genome_start = node_end->position;
    int lim_ref;
    lim_ref = genome_end - genome_start;	    
    if (lim_ref > gap_close) {
      lim_ref = gap_close;
      genome_start = genome_end - lim_ref - 1;
    }

    //printf("lim_ref(%i) = genome_end(%lu) - genome_start(%lu) | gap_close = %i\n", lim_ref, genome_end, genome_start, gap_close);
    if (lim_ref < 0) {
      //For first step, genome_start = first_cal->end
      int dsp = abs(lim_ref);
      gap_close += dsp;
      read_pos  += dsp;
      array_list_insert((void *)genome_start + 1, final_positions);
      //printf(" Insert 1) %lu\n", genome_start + 1);
      //printf(":::--::::GENOME GAP %i\n", genome_end);
      //printf("RECALCULATING GAP AND READ POS (%i) (gap_close, read_pos)(%i, %i)\n", dsp, gap_close, read_pos);
    } else {
      assert((genome_end - genome_start) + 2 < 2048);

      genome_read_sequence_by_chr_index(reference, 0, 
					cal_chromosome_id - 1, 
					&genome_start, &genome_end, genome);
      //printf("(%lu-%lu)Reference START_SP [[[START_SP]]]----[END_SP]: %s\n",  genome_start, genome_end, reference);
      int t, c;
      //printf("Num targets %i\n", num_targets);
      for (t = 0; t < num_targets; t++) {
	node_end = array_list_get(t, ends_targets);
	lim_ref = genome_end - node_end->position;

	if (lim_ref > gap_close) { 
	  lim_ref = gap_close;
	}
	dist = 0;
	//printf("\tTravel to %lu\n", node_end->position);
	reference_len = strlen(reference) - 1;
	for (c  = 0; c < lim_ref; c++) { 
	  //printf("\t\t[%c vs %c]\n", query_map[read_pos - c], reference[reference_len - c]);
	  if ((read_pos - c) < 0 || (reference_len - c) < 0) { goto sp_free; }
	  if (query_map[read_pos - c] != reference[reference_len - c]) { 
	    dist++;
	  }
	}
	//max_dist = lim_ref <= 5 ? 2 : lim_ref / 4 + 1;
	max_dist = lim_ref <= 5 ? 2 : lim_ref / 4;
	//printf("dist = %i < max_dist = %i\n", dist, max_dist);
	if (dist <= max_dist) { 
	  //printf("\t\tMAX DISTANCE GOOD %i!\n", dist);
	  //max_dist += dist;
	  final_dist += dist;
	  if (lim_ref == gap_close) { t++; break; }
	} else {
	  break;
	}
      }
      t--;
      if (t < 0) { break; }
      else {
	node_end = array_list_get(t, ends_targets);
	lim_ref = genome_end - node_end->position;
	if (lim_ref > gap_close) { lim_ref = gap_close; }
	gap_close -= lim_ref;
	read_pos -= lim_ref;
	//printf("GAP CLOSE (%i) READ POS (%i)\n", gap_close, read_pos);
	if (gap_close <= 0) {
	  //Final Map
	  map = 1;
	  array_list_insert((void *)genome_end - c + 1, final_positions);
	  //printf(" Insert 2) %lu\n", genome_end - c + 1);
	  break;
	} else {
	  array_list_insert((void *)node_end->position + 1, final_positions);
	  //printf(" Insert 3) %lu\n", node_end->position + 1);
	}
      }
    }
	    
    //3.Select the correct start
    int pos;
    array_list_t *starts_list = ((end_data_t *)node_end->data)->list_starts;

    array_list_clear(final_ends, (void *)NULL);
    //Filter ends before cal_next->start
    for (int e = 0; e < array_list_size(starts_list); e++) {
      size_t start = (size_t)array_list_get(e, starts_list);
      if (filter_pos) {
	if (start > filter_pos) {
	  array_list_insert((void *)start, final_ends);
	}
      } else {
	array_list_insert((void *)start, final_ends);
      }
    }
    if (array_list_size(final_ends) >= 1) {
      //Select best start
      char s_reference[array_list_size(final_ends)][2048];
      int references_len[array_list_size(final_ends)];
      //Making reference...
      for (int s = 0; s < array_list_size(final_ends); s++) {
	genome_end = (size_t)array_list_get(s, final_ends) - 1;
	genome_start = genome_end - gap_close - 1;
	assert((genome_end - genome_start) + 2 < 2048);
	genome_read_sequence_by_chr_index(s_reference[s], 0, 
					  cal_chromosome_id - 1, 
					  &genome_start, &genome_end, genome); 	      
	references_len[s] = strlen(s_reference[s]) - 1;
	//printf("(%lu-%lu)Reference END_SP [START_SP]----[[[END_SP]]](%i): %s\n", genome_start, genome_end, s, s_reference[s]);
      } 
      //Select the best start
      //int lim = (array_list_get(1, final_starts)) - (array_list_get(0, final_starts));
      int const MIN_ANCHOR = 20;
      int anchor_nt = MIN_ANCHOR; 
      if (anchor_nt > gap_close) { 
	anchor_nt = gap_close;
      }
      pos = 0;
      int *ref_matches = (int *)calloc(array_list_size(final_ends), sizeof(int));
      int *ref_mismatches = (int *)calloc(array_list_size(final_ends), sizeof(int));
      //printf("anchor nt = %i, array_list_size(final_ends) = %i\n", anchor_nt, array_list_size(final_ends));
      for (int c = 0; c < anchor_nt; c++) {
	for (int s = 0; s < array_list_size(final_ends); s++) {	  
	  if (read_pos - c < 0 || references_len[s] - c < 0) { goto sp_free; }
	  //printf("\t[%c vs %c]\n", query_map[read_pos - c], s_reference[s][references_len[s] - c]);
	  if (query_map[read_pos - c] != s_reference[s][references_len[s] - c]) {
	    ref_mismatches[s]++;
	  } else {
	    ref_matches[s]++;
	  }
	}
      }
      float score, max_score = 0.0;
      for (int s = 0; s < array_list_size(final_ends); s++) {
	score = ref_matches[s]*0.5 - ref_mismatches[s]*0.4;
	if (score > max_score) { 
	  max_score = score;
	  pos = s;
	}
      }		
      max_dist = 3;
      if (ref_mismatches[pos] > max_dist) {
	free(ref_matches);
	free(ref_mismatches);
	break;
      }

      final_dist += ref_mismatches[pos]; 

      final_pos = (size_t)array_list_get(pos, final_ends);
      array_list_insert((void *)final_pos, final_positions);
      //printf(" Insert 4) %lu\n", final_pos);
      genome_end = final_pos - anchor_nt - 1;
      gap_close -= anchor_nt;
      read_pos -= anchor_nt;
      //printf("FINAL MATCHES (%i) VS FINAL MISMATCHES (%i)\n", ref_matches[s], ref_mismatches[s]);
      //printf("GAP CLOSE (%i) READ POS (%i)\n", gap_close, read_pos);

      free(ref_mismatches);
      free(ref_matches);

      if (gap_close <= 0) {
	array_list_insert((void *)final_pos - anchor_nt, final_positions);
	//printf(" Insert 5) %lu\n", final_pos - anchor_nt);
	map = 1;
	break;
      } else {
	//printf("Search in meta: %lu-%lu\n", final_pos, final_pos - anchor_nt);
	if (metaexon_search(cal_strand, cal_chromosome_id - 1, 
			    final_pos - anchor_nt, final_pos, 
			    &final_metaexon, metaexons)) {
	  if (final_metaexon) {
	    if (final_metaexon->left_closed) {
	      array_list_clear(ends_targets, (void *)NULL);
	      for (int s_0 = 0; s_0 < array_list_size(final_metaexon->left_breaks); s_0++) {
		node_end = array_list_get(s_0, final_metaexon->left_breaks);
		if (final_pos + 10 >= node_end->position) {
		  //printf("Insert node_end\n");
		  array_list_insert(node_end, ends_targets);		    
		}
	      }
	    } else {
	      //We can close the gap
	      int exon_length = final_metaexon->end - final_metaexon->start;
	      dist = 0;
	      if (exon_length >= gap_close) {
		genome_start = genome_end - gap_close - 1;
		assert((genome_end - genome_start) + 2 < 2048);
		assert((read_pos) <= strlen(query_map));

		genome_read_sequence_by_chr_index(reference, 0, 
						  cal_chromosome_id - 1, 
						  &genome_start, &genome_end, genome);		
		//printf("CLOSE GAP (%i): %s\n", gap_close, reference);
		reference_len = strlen(reference) - 1;
		for (int c  = 0; c < gap_close; c++) { 
		  //printf("\t[%c vs %c]\n", query_map[read_pos - c], reference[reference_len - c]);
		  if ((read_pos - c < 0) &&
		      (reference_len - c < 0)) { goto sp_free; }

		  if (query_map[read_pos - c] != reference[reference_len - c]) { 
		    dist++;
		  }
		}
		if (dist < 5) {
		  array_list_insert((void *)final_pos - gap_close - anchor_nt, final_positions);
		  //printf(" Insert 6) %lu\n", final_pos - gap_close - anchor_nt);
		  final_dist += dist; 
		  map = 1;
		}
	      }
	      //printf("We close gap! with  map=%i\n", map);
	      break;
	    }
	  }
	} else {
	  //printf("NOT  FOUND\n");
	  break;
	}
      }
    } else {
      //printf("Not final ends :(\n");
      break;
    }
  }

  if (map) {
    size_t pos_prev = cal->end, pos_next;
    cigar_op_t *op;
    cigar_code = cigar_code_new();
    cigar_code->distance = final_dist;

    pos_prev = cal->start - 1;
    for (int sp = 0; sp < array_list_size(final_positions); sp++) {
      pos_next = (size_t)array_list_get(sp, final_positions);
      //printf("%lu - %lu(%i/%i)\n", pos_prev, pos_next, sp, array_list_size(final_positions));
      if (sp % 2 == 0) {		
	op = cigar_op_new(pos_prev - pos_next + 1, 'M');
      } else {
	//int aux = pos_prev - pos_next + 1;
	//op = cigar_op_new(aux, 'N');

	size_t start_sp = pos_next;
	size_t end_sp = pos_prev;
	int aux = end_sp - start_sp + 1;

	op = cigar_op_new(aux, 'N');

	char nt_end[5], nt_start[5];
	int type;
	//printf("%lu - %lu(%i/%i) = %i\n", pos_prev, pos_next, sp, array_list_size(final_positions), aux);	
	//Report splice junction
	size_t g_start = start_sp;
	size_t g_end   = start_sp + 1;	
	genome_read_sequence_by_chr_index(nt_start, cal->strand, cal->chromosome_id - 1, &g_start, &g_end, genome);

	g_end   = end_sp;
	g_start = end_sp - 1;
	genome_read_sequence_by_chr_index(nt_end, cal->strand, cal->chromosome_id - 1, &g_start, &g_end, genome);
		
	type = splice_junction_type(nt_start[0], nt_start[1], nt_end[0], nt_end[1]);
	int splice_strand;
	avl_node_t *avl_node_start, *avl_node_end;
	
	splice_strand = 0;
	if (type == CT_AC_SPLICE || type == GT_AT_SPLICE || type == CT_GC_SPLICE ) {
	  splice_strand = 1;
	}

	//printf("R:[%i:%lu-%lu] : %c%c vs %c%c\n", cal->chromosome_id, start_sp, end_sp, nt_start[0], nt_start[1], nt_end[0], nt_end[1]);

	char str_sp_type[10];
	if (type == NOT_SPLICE) {
	  type = UNKNOWN_SPLICE;
	  sprintf(str_sp_type, "%c%c-%c%c\0", nt_start[0], nt_start[1], nt_end[0], nt_end[1]);
	}

	
	allocate_start_node(cal->chromosome_id - 1,
			    splice_strand,
			    start_sp,
			    end_sp,
			    start_sp,
			    end_sp,
			    FROM_READ,
			    type,
			    str_sp_type, 
			    &avl_node_start,
			    &avl_node_end, 
			    avls_list);
      }
      //printf("\tADD %i%c\n", op->number, op->name);
      cigar_code_insert_first_op(op, cigar_code);
      pos_prev = pos_next - 1;
    }

    //printf(" :::::::::::::::::::::::::::: %s\n", new_cigar_code_string(cigar_code));
    //cal->info = cigar_code;
    //s_prev->info = cigar_code;
    //meta_alignment = meta_alignment_new();
    //meta_alignment_insert_cal(cal, meta_alignment);
    //meta_alignment_set_status(META_CLOSE, meta_alignment);
	      
  }

 sp_free:
  array_list_free(final_positions, NULL);
  array_list_free(ends_targets, NULL);
  array_list_free(final_ends, NULL);

  //if (cigar_code == NULL) { exit(-1); }

  //printf("Final Cigar %s\n", new_cigar_code_string(cigar_code));

  return cigar_code;

}

cigar_code_t *search_double_anchors_cal(char *query_map,
					cal_t *first_cal, 
					cal_t *last_cal,
					metaexons_t *metaexons, 
					genome_t *genome,
					fastq_read_t *fq_read,
					int *type, 
					avls_list_t *avls_list) {
  //return NULL;
  int l_found = 0;
  int r_found = 0;
  int read_nt, read_orig_nt;
  metaexon_t *first_metaexon = NULL, *last_metaexon = NULL;
  assert(first_cal != NULL);
  assert(last_cal != NULL);

  seed_region_t *s_prev = linked_list_get_last(first_cal->sr_list);    
  seed_region_t *s_next = linked_list_get_first(last_cal->sr_list);	
  assert(s_prev != NULL);
  assert(s_next != NULL);

  //int found_sp = 0;
  size_t genome_start, genome_end;
  char reference[2048];
  //meta_alignment_t *meta_alignment = NULL;
  cigar_code_t *cigar_code = NULL;

  *type = META_ALIGNMENT_MIDDLE;

  pthread_mutex_lock(&metaexons->mutex[first_cal->chromosome_id - 1]);

  //printf("-------------------------------------------\n");
  //cal_print(first_cal);
  //cal_print(last_cal);
  //printf("...........................................\n");
  //metaexons_print(metaexons);
  //printf("-------------------------------------------\n");

  int m_found_f = metaexon_search(first_cal->strand, 
				  first_cal->chromosome_id - 1,
				  first_cal->start,
				  first_cal->end, &first_metaexon,
				  metaexons); 

  int m_found_l = metaexon_search(last_cal->strand,
				  last_cal->chromosome_id - 1,
				  last_cal->start, 
				  last_cal->end,
				  &last_metaexon,
				  metaexons);

  //if (first_metaexon == last_metaexon) {
  //return NULL;
  //}

  if (m_found_f) {
    //printf("First Meta %i\n", first_metaexon->right_closed);
    l_found = first_metaexon->right_closed;
  } //else {
    //printf("Not found\n");
  //}	 

  if (m_found_l) {
    //printf("Last Meta %i\n", last_metaexon->left_closed);
    r_found = last_metaexon->left_closed;
    //} else {
    //printf("Not found\n");
  }
  
  //printf("l_found = %i, r_found = %i\n", l_found, r_found);
  if (l_found && r_found) {
    if (last_metaexon->start < first_metaexon->start) {
      cal_print(first_cal);
      cal_print(last_cal);
      LOG_FATAL_F("first_meta(%p): [%lu-%lu], last_meta(%p): [%lu - %lu]\n", 
		  first_metaexon,
		  first_metaexon->start, 
		  first_metaexon->end,
		  last_metaexon,
		  last_metaexon->start, 
		  last_metaexon->end);
    }
    
    //Is the correct splice?? Ask to length exons...
    //printf("Found Metaexon.... Is the correct???\n");	    
    array_list_t *introns_targets = array_list_new(10, 1.5f, COLLECTION_MODE_ASYNCHRONIZED);
    array_list_t *starts_targets  = array_list_new(10, 1.5f, COLLECTION_MODE_ASYNCHRONIZED);
    array_list_t *ends_targets  = array_list_new(10, 1.5f, COLLECTION_MODE_ASYNCHRONIZED);
    intron_t *intron;
    int size_ex_l, size_ex_r, dsp_l, dsp_r;
    avl_node_t *node_start, *node_end;

    //1st Search posible starts splice junctions
    //if (array_list_size(first_metaexon->right_breaks) > 200) {
      /*    printf("OVERFLOW 1 [%lu-%lu]:%i\n", first_metaexon->start, first_metaexon->end, 
	     array_list_size(first_metaexon->right_breaks));

      if (metaexons->bypass_pointer[first_cal->chromosome_id - 1][first_cal->start / 1024].first) {
	linked_list_item_t *list_item = metaexons->bypass_pointer[first_cal->chromosome_id - 1][first_cal->start / 1024].first;
	int n = 0, nMax = 10;
	if (list_item) {
	  metaexon_t *metaexon = list_item->item;
	  while (n < nMax) {
	    printf("[%lu - %lu]-", metaexon->start, metaexon->end);
	    list_item = list_item->next;
	    if (list_item != NULL) {
	      metaexon = list_item->item;
	    } else { 
	      break;
	    }
	    n++;
	  }
	  printf("\n");
	}
      }
      */
      //return NULL;
      //exit(-1);
      //}

    for (int s_0 = 0; s_0 < array_list_size(first_metaexon->right_breaks); s_0++) {
      node_start = array_list_get(s_0, first_metaexon->right_breaks);
      //printf("\t\t Node start %lu\n", node_start->position);
      if (node_start == NULL) { printf("node start NULL\n"); exit(-1); }
      if (first_cal->end - 10 <= node_start->position) {
	array_list_insert(node_start, starts_targets);
      }
    }
    
    //if (array_list_size(last_metaexon->left_breaks) > 200) {
      /*printf("OVERFLOW 2 [%lu-%lu]:%i\n", last_metaexon->start, last_metaexon->end, 
	     array_list_size(last_metaexon->left_breaks));

      if (metaexons->bypass_pointer[last_cal->chromosome_id - 1][last_cal->start / 1024].first) {
	linked_list_item_t *list_item = metaexons->bypass_pointer[last_cal->chromosome_id - 1][last_cal->start / 1024].first;
	int n = 0, nMax = 10;
	if (list_item) {
	  metaexon_t *metaexon = list_item->item;
	  while (n < nMax) {
	    printf("[%lu - %lu]-", metaexon->start, metaexon->end);
	    list_item = list_item->next;
	    if (list_item != NULL) {
	      metaexon = list_item->item;
	    } else { 
	      break;
	    }
	    n++;
	  }
	  printf("\n");
	}
      }
      */
      //return NULL;
      //exit(-1);
    //}

    //assert(last_metaexon->left_);
    for (int e_0 = 0; e_0 < array_list_size(last_metaexon->left_breaks); e_0++) {
      node_end = array_list_get(e_0, last_metaexon->left_breaks);
      if (node_end == NULL) { printf("node end NULL\n"); exit(-1); }
      //printf("\t\t Node end %lu\n", node_end->position);
      if (last_cal->start + 10 >= node_end->position) {
	array_list_insert(node_end, ends_targets);
      }
    } 
    
    for (int s_0 = 0; s_0 < array_list_size(starts_targets); s_0++) {
      node_start = array_list_get(s_0, starts_targets);
      array_list_t *ends_list = ((start_data_t *)node_start->data)->list_ends;
      for (int e_0 = 0; e_0 < array_list_size(ends_list); e_0++) {
	splice_end_t *splice_end = array_list_get(e_0, ends_list);
	for (int e_1 = 0; e_1 < array_list_size(ends_targets); e_1++) {
	  avl_node_t *node_aux = array_list_get(e_1, ends_targets);
	  //printf("\t\t Node end %lu == %lu\n", splice_end->end, node_aux->position);
	  if (splice_end->end == node_aux->position) {
	    //*******************************************//
	    //     [  L_EXON.1  ]-------[  R_EXON.1  ]   //
	    //     [  L_EXON.2 ]----------[R_EXON.2  ]   //
	    //     [  L_EXON.3]-------------[R_EXON.3]   //
	    //*******************************************//
	    if (node_start->position > splice_end->end) {
	      LOG_FATAL_F("%lu - %lu intron ERROR!\n", node_start->position, splice_end->end);
	    }
	    intron = intron_new(last_cal->strand, last_cal->chromosome_id, 
				node_start->position, splice_end->end);
	    array_list_insert(intron, introns_targets);
	  }
	}
      }
    }
	  
    read_nt = s_next->read_start - s_prev->read_end - 1; //ok
    //if (read_nt < 0) { read_nt = 0; }	    

    read_orig_nt = read_nt;
    //We not found introns... We have one exon between CALs??  //min-exon 30nt
    
    if (array_list_size(introns_targets) <= 0) { 
      //printf("INTRON TARGETS, SEARCH LEFT\n");
      read_nt = fq_read->length - (s_prev->read_end + 1);
      cigar_code = search_left_single_anchor(read_nt, 
					     first_cal,
					     0,
					     first_metaexon->right_breaks,
					     query_map, metaexons, genome, 
					     avls_list);
      *type = META_ALIGNMENT_LEFT;
    }    

    //metaexons_show(metaexons);
    //assert(array_list_size(introns_targets) != 0);
    //printf("read_gap = (%i - %i) %i, introns_targets = %i\n", 
    //	   s_next->read_start, s_prev->read_end, read_nt, array_list_size(introns_targets));
    //return NULL;

    size_t s_intron, e_intron;
    int distance = 0;
    size_t s_prev_read_end, s_prev_genome_end, first_cal_end;
    size_t s_next_read_start, s_next_genome_start, last_cal_start;

    //printf("read_nt = %i\n", read_nt);
    //Check introns... good luck! :)
    for (int in = 0; in < array_list_size(introns_targets); in++) {
      intron = array_list_get(in, introns_targets);
      assert(intron != NULL);
      //printf("Intron target %i:[%lu-%lu] vs cal[%lu-%lu]\n", in, intron->start, intron->end,
      //     first_cal->end, last_cal->start);
      
      read_nt             = read_orig_nt;
      s_prev_read_end     = s_prev->read_end;
      s_prev_genome_end   = s_prev->genome_end;
      first_cal_end       = first_cal->end;
      s_next_read_start   = s_next->read_start; 
      s_next_genome_start = s_next->genome_start; 
      last_cal_start      = last_cal->start; 

      dsp_l = 0;
      dsp_r = 0;

      size_ex_l = (int)(intron->start - (first_cal->end + 1));
      //printf("\t LEFT = (i)%lu - %lu = %i\n", intron->start, first_cal->end + 1, size_ex_l);
      if (size_ex_l < 0) {
	dsp_l = abs(size_ex_l);
	s_prev_read_end -= dsp_l;
	s_prev_genome_end -= dsp_l;
	first_cal_end -= dsp_l;
	read_nt += dsp_l;
	size_ex_l = 0;
      }
	    
      size_ex_r = (int)(last_cal->start - (intron->end + 1));
      //printf("\t RIGHT = %lu - (i)%lu = %i\n", last_cal->start, intron->end + 1, size_ex_r);
      if (size_ex_r < 0) { 
	dsp_r = abs(size_ex_r); 
	s_next_read_start += dsp_r;
	s_next_genome_start += dsp_r;
	last_cal_start += dsp_r;
	read_nt += dsp_r;
	size_ex_r = 0;
      }
      
      if (read_nt < 0) { printf("ERROR: read_nt = %i (%s)\n", read_nt, fq_read->id); cal_print(first_cal); cal_print(last_cal); exit(-1); }

      //assert(read_nt >= 0);
      
      //printf("size_ex_l = %i, size_ex_r = %i, read_nt = %i\n", size_ex_l, size_ex_r, read_nt);
      if (read_nt == (size_ex_l + size_ex_r)) {
	//Close gap		
	if (size_ex_l > 1) {
	  genome_start = first_cal_end + 1;
	  genome_end   = first_cal_end + size_ex_l;
	  
	  if ((genome_end - genome_start + 5) >= 2048) {
	    LOG_FATAL_F("%lu - %lu\n", genome_start, genome_end);
	  }	  

	  assert((s_prev_read_end + size_ex_l + 1) < strlen(query_map));

	  genome_read_sequence_by_chr_index(reference, 0, 
					    first_cal->chromosome_id - 1, 
					    &genome_start, &genome_end, genome);
	  //printf("Ref-l-%s\n", reference);
	  for (int c = 0; c < size_ex_l; c++) {
	    //printf("%c/%c | ", query_map[s_prev->read_end + 1 + c], reference[c]);
	    //if ((s_prev_read_end + 1 + c) >= strlen(query_map)) { exit(-1); }

	    if (query_map[s_prev_read_end + 1 + c] != reference[c]) {
	      distance++;
	    }
	  }
	  //printf("\n");
	} else if (size_ex_l == 1) {
	  distance++;
	}
	//printf("1st-Distance-l %i\n", distance);
	//return NULL;

	s_prev_read_end += size_ex_l;
	s_prev_genome_end += size_ex_l;
	first_cal_end += size_ex_l;
	
	if (size_ex_r > 1) {
	  genome_start = last_cal_start - size_ex_r;
	  genome_end   = last_cal_start;
	  if (genome_end - genome_start + 5 >= 2048) {
	    printf("%lu - %lu = %lu\n", genome_start, genome_end, genome_end - genome_start);
	    exit(-1);
	  }

	  genome_read_sequence_by_chr_index(reference, 0, 
					    first_cal->chromosome_id - 1, 
					    &genome_start, &genome_end, genome);		  
	  //printf("Ref-r-%s\n", reference);
	  for (int c = 0; c < size_ex_r; c++) {
	    //printf("%c/%c | ", query_map[s_next->read_start - size_ex_r + c], reference[c]);
	    if ((s_next_read_start - size_ex_r + c) >= strlen(query_map) || 
		(s_next_read_start - size_ex_r + c) < 0) {
	      printf("ERROR DOUBLE ANCHOR OVERFLOW\n");
	      exit(-1);
	    }
	    if (query_map[s_next_read_start - size_ex_r + c] != reference[c]) {
	      distance++;
	    }
	  }
	  //printf("\n");
	}else if (size_ex_r == 1) {
	  distance++;
	}

	s_next_read_start -= size_ex_r;
	s_next_genome_start -= size_ex_r;
	last_cal_start -= size_ex_r;
	
	//seed_region = seed_region_new(s_next->read_start, s_next->read_end, 
	//			      s_next->genome_start, s_next->genome_end, 1);
	//linked_list_insert_last(seed_region, first_cal->sr_list);		
	//first_cal->end = last_cal->end;
	//printf("2nd-Distance-r %i\n", distance);
	if (distance > (size_ex_l + size_ex_r)/2) {
	  continue;
	}

	cigar_code = cigar_code_new();

	size_ex_l -= dsp_l;
	size_ex_r -= dsp_r;

	cigar_code_append_new_op(size_ex_l, 'M', cigar_code);
	cigar_code_append_new_op(last_cal_start - first_cal_end - 1, 'N', cigar_code);
	cigar_code_append_new_op(size_ex_r, 'M', cigar_code);		
	cigar_code->distance = distance;

	//break;
	//Allocate splice
	size_t start_sp = first_cal_end + 1;
	size_t end_sp = last_cal_start - 1;
	int aux = end_sp - start_sp;
	char nt_end[10], nt_start[10];
	int type;

	//Report splice junction
	size_t g_start = start_sp;
	size_t g_end   = start_sp + 1;	
	genome_read_sequence_by_chr_index(nt_start, first_cal->strand, first_cal->chromosome_id - 1, &g_start, &g_end, genome);

	g_end   = end_sp;
	g_start = end_sp - 1;
	genome_read_sequence_by_chr_index(nt_end, first_cal->strand, first_cal->chromosome_id - 1, &g_start, &g_end, genome);
		
	type = splice_junction_type(nt_start[0], nt_start[1], nt_end[0], nt_end[1]);
	int splice_strand;
	avl_node_t *avl_node_start, *avl_node_end;
	
	splice_strand = 0;
	if (type == CT_AC_SPLICE || type == GT_AT_SPLICE || type == CT_GC_SPLICE ) {
	  splice_strand = 1;
	}

	//printf("M:[%i:%lu-%lu] : %c%c vs %c%c\n", first_cal->chromosome_id, start_sp, end_sp, nt_start[0], nt_start[1], nt_end[0], nt_end[1]);

	char str_sp_type[10];
	if (type == NOT_SPLICE) {
	  type = UNKNOWN_SPLICE;
	  sprintf(str_sp_type, "%c%c-%c%c\0", nt_start[0], nt_start[1], nt_end[0], nt_end[1]);
	}

	allocate_start_node(first_cal->chromosome_id - 1,
			    splice_strand,
			    start_sp,
			    end_sp,
			    start_sp,
			    end_sp,
			    FROM_READ,
			    type,
			    str_sp_type, 
			    &avl_node_start,
			    &avl_node_end, 
			    avls_list);

	break; //If found break loop	
      } else {
	LOG_DEBUG_F("@@@PROBLEM CALCULATING SPLICE JUNCTION id = %s\n", fq_read->id);
      }		
      //printf("Distance to close start splice %i\n", size_ex_l);
      //printf("Distance to close end   splice %i\n", size_ex_r);	      
    } //End loop introns	    
    array_list_free(starts_targets, (void *)NULL);
    array_list_free(ends_targets, (void *)NULL);
    array_list_free(introns_targets, (void *)intron_free);
  }//End if (l_found && r_found)

  //if (cigar_code == NULL) { exit(-1); }

  pthread_mutex_unlock(&metaexons->mutex[first_cal->chromosome_id - 1]);  

  return cigar_code;

}


//IF MODE == 0, DOUBLE ANCHORS
//IF MODE == 1, LEFT ANCHOR
//IF MODE == 2, RIGHT ANCHOR
int generate_cals_between_anchors (int mode,
				   int num_chromosomes,
				   cal_t *first_cal, 
				   cal_t *last_cal, 
				   fastq_read_t *fq_read, 
				   array_list_t *seeds_list, 
				   array_list_t *cals_list, 
				   cal_optarg_t *cal_optarg) {
  //printf("GENERATE CALS\n");
  seed_region_t *s_prev, *s_next;
  int min_seeds = 0, max_seeds = 1000;
  region_t *region_prev = NULL, *region_next = NULL;
  array_list_t *mapping_list = array_list_new(500, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  int strand, chromosome_id;
  int seed_size = cal_optarg->seed_size;
  const int max_intron_size = cal_optarg->max_intron_size;

  array_list_clear(cals_list, (void *)NULL);

  if (mode == CALING_DOUBLE_ANCHORS || 
      mode == CALING_LEFT_ANCHORS) {
    s_prev = linked_list_get_last(first_cal->sr_list);    
    region_prev = region_bwt_new(first_cal->chromosome_id,
				 first_cal->strand,
				 first_cal->start, 
				 first_cal->end,
				 s_prev->read_start,
				 s_prev->read_end,
				 fq_read->length,
				 0);
    array_list_insert(region_prev, mapping_list);
    strand = first_cal->strand;
    chromosome_id = first_cal->chromosome_id;
  }

  //printf("Region FORW(%i)[%i:%lu|%i-%i|%lu]\n ", region_prev->id, region_prev->chromosome_id, 
  //	   region_prev->start, region_prev->seq_start, region_prev->seq_end, region_prev->end);	  
  if (mode == CALING_DOUBLE_ANCHORS || 
      mode == CALING_RIGHT_ANCHORS) {    
    s_next = linked_list_get_first(last_cal->sr_list);
    int id, gap;
    if (mode == CALING_DOUBLE_ANCHORS) {
      gap = s_next->read_start - s_prev->read_end;
      id = (gap / seed_size) + ((gap % seed_size) > 0) + 1;
    } else {
      id = 1;
    }
    region_next = region_bwt_new(last_cal->chromosome_id,
				 last_cal->strand,
				 last_cal->start, 
				 last_cal->end,
				 s_next->read_start,
				 s_next->read_end,
				 fq_read->length,
				 id);
    array_list_insert(region_next, mapping_list);	      
    strand = last_cal->strand;
    chromosome_id = last_cal->chromosome_id;
  }
  //printf("Region BACK(%i)[%i:%lu|%i-%i|%lu]\n ", region_next->id, region_next->chromosome_id, 
  //	   region_next->start, region_next->seq_start, region_next->seq_end, region_next->end);  

  //printf(":::::::::::::: seeds %i\n", array_list_size(seeds_list));
  for (int j = 0; j < array_list_size(seeds_list); j++) {
    region_t *region = array_list_get(j, seeds_list);

    //printf("@@@@@@@(%i)Region [%i:%lu|%i-%i|%lu]: ", region->id, region->chromosome_id, 
    //	     region->start, region->seq_start, region->seq_end, region->end);

    if (region->strand == strand && 
	region->chromosome_id == chromosome_id) {
      if (mode == CALING_DOUBLE_ANCHORS) {
	if (region_prev->end - seed_size <= region->start &&
	    region_next->start + seed_size >= region->end) {
	  //printf("Insert\n");
	    array_list_insert(region, mapping_list);
	} else { /*printf("Not Insert\n");*/ }
      } else if (mode == CALING_LEFT_ANCHORS) {
	if (region_prev->end - seed_size <= region->start) {
	  //printf("Insert\n");
	  array_list_insert(region, mapping_list);
	} else { /*printf("Not Insert\n");*/ } 
      } else {
	if (region_next->start + seed_size >= region->end) {
	  //printf("Insert\n");
	  array_list_insert(region, mapping_list);
	} else { /*printf("Not Insert\n");*/ }
      }
    } else {
      //printf("Not Insert\n");
    }
  }

    
  int min_cal_size = cal_optarg->min_cal_size;
  if (mode == CALING_DOUBLE_ANCHORS ) {
    if (last_cal->start - first_cal->end < max_intron_size) {
      min_cal_size = 15;      
    }
  }
  

  bwt_generate_cal_list_linked_list(mapping_list,
				    cal_optarg,
				    &min_seeds, &max_seeds,
				    num_chromosomes,
				    cals_list,
				    fq_read->length,
				    min_cal_size, 0);

  //cal_optarg->min_cal_size = min_cal_size;
  int num_cals = array_list_size(cals_list);
  if (num_cals > 0) {
    for (int j = num_cals - 1; j >= 0; j--) {
      cal_t *cal = array_list_get(j, cals_list);
      seed_region_t *s = linked_list_get_first(cal->sr_list);
      if (s == NULL) {
	array_list_remove_at(j, cals_list);
	cal_free(cal);
      } else {
	if (s->genome_start != cal->start) { 
	  array_list_remove_at(j, cals_list);
	  cal_free(cal);
	} else {
	  seed_region_t *s = linked_list_get_last(cal->sr_list);
	  if (s == NULL) {
	    array_list_remove_at(j, cals_list);
	    cal_free(cal);
	  } else {
	    if (s != NULL && s->genome_end != cal->end) { 
	      array_list_remove_at(j, cals_list);
	      cal_free(cal);
	    }
	  }
	}
      }
    }
  }

  num_cals = array_list_size(cals_list);

  int founds[num_cals];
  int found = 0;
  cal_t *cal;
  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cals_list);
    LOG_DEBUG_F("\tcal %i of %i: sr_list size = %i (cal->num_seeds = %i) %i:%lu-%lu\n", 
		j, num_cals, cal->sr_list->size, cal->num_seeds,
		cal->chromosome_id, cal->start, cal->end);
    founds[j] = 0;
    if (cal->sr_list->size > 0) {
      int start = 0;
      for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
	seed_region_t *s = list_item->item;		  
	LOG_DEBUG_F("\t\t:: star %lu > %lu s->read_start\n", start, s->read_start);
	if (start > s->read_start || 
	    s->read_start >= s->read_end) {
	  LOG_DEBUG("\t\t\t:: remove\n");
	  found++;
	  founds[j] = 1;
	}
	start = s->read_end + 1;
      }
    } else {
      found++;
      founds[j] = 1;
    }
  }

  //Check CALs TODO:Delete?Duplicate??
  int start = 0;
  for (int c = 0; c < num_cals; c++) { 
    cal = array_list_get(c, cals_list);
    seed_region_t *s_prev_aux = linked_list_get_first(cal->sr_list);    
    seed_region_t *s_next_aux = linked_list_get_last(cal->sr_list);
    if (s_prev_aux == NULL || s_next_aux == NULL) { 
      founds[c] = 1;
      found++;
    } else {
      //printf("CAL %i: %i < %i\n", c, s_prev_aux->read_start, start);
      if (s_prev_aux->read_start < start) {
	founds[c] = 1;
	found++;	
      } else {
	start = s_next_aux->read_end;
      }
    }
  }
	      
  if (found) {
    for (int j = num_cals - 1; j >= 0; j--) {
      if (founds[j]) {
	cal = array_list_remove_at(j, cals_list);
	array_list_free(cal->candidates_seeds_start, NULL);
	array_list_free(cal->candidates_seeds_end, NULL);
	if (cal->sr_list != NULL) { 
	  linked_list_free(cal->sr_list, NULL); 
	}
	if (cal->sr_duplicate_list != NULL) { 
	  linked_list_free(cal->sr_duplicate_list, NULL); 
	}
	free(cal);
	//array_list_set(j, NULL, cals_list);
	//cal_free(cal);
      }
    }
  }
  //seeds 
  
  cal_free(first_cal);
  cal_free(last_cal);

  //printf(" ==== FINAL CALS (%i) ==== \n", array_list_size(cals_list));
  //for (size_t j = 0; j < array_list_size(cals_list); j++) {
  //cal = array_list_get(j, cals_list);
  //cal_print(cal);
  //}
  //printf(" ==== END FINAL CALS ==== \n");
  array_list_free(mapping_list, NULL);

  if (region_prev) { region_bwt_free(region_prev); }
  if (region_next) { region_bwt_free(region_next); }

  return array_list_size(cals_list);
}

int meta_alignment_regenerate_cigar(meta_alignment_t *meta_alignment, char *query_map, genome_t *genome) {
  array_list_t *cals_list = meta_alignment->cals_list;

  /*
  printf("==== CALS LIST %i ====\n", array_list_size(cals_list));
  for (int j = 0; j < array_list_size(cals_list); j++) {
    cal_t *cal = array_list_get(j, cals_list);
    cal_print(cal);
  }
  printf("==== CALS LIST ====\n");      
  */

  int num_cals = array_list_size(cals_list);
  if (!num_cals) { return 1; }

  cigar_code_t *cigar_code_n = cigar_code_new();
  cal_t *cal_prev = array_list_get(0, cals_list);
  cigar_code_t *cigar_code_prev = cal_prev->info;

  if (meta_alignment->cigar_left != NULL) {
    cigar_code_t *cigar_code_aux = meta_alignment->cigar_left;
    for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
      cigar_op_t *op = array_list_get(i, cigar_code_aux->ops);
      cigar_code_append_new_op(op->number, op->name, cigar_code_n);
    } 
    cigar_code_n->distance += cigar_code_aux->distance;
  }

  for (int i = 0; i < array_list_size(cigar_code_prev->ops); i++) {
    cigar_op_t *op = array_list_get(i, cigar_code_prev->ops);
    cigar_code_append_new_op(op->number, op->name, cigar_code_n);
  }

  cal_t *cal_next;
  int read_gap;
  seed_region_t *seed_prev, *seed_next;
  int read_start, read_end;
  size_t genome_start, genome_end;
  char reference_prev[2048];
  char reference_next[2048];
  int max_distance = 3;
  int distance;
  int prev_M, next_M;

  if (num_cals > 1) {
    for (int i = 1; i < num_cals; i++) {
      cal_next = array_list_get(i, cals_list);
      seed_prev = linked_list_get_last(cal_prev->sr_list);
      seed_next = linked_list_get_first(cal_next->sr_list);

      read_gap = seed_next->read_start - seed_prev->read_end - 1;

      read_start = seed_prev->read_end + 1;
      read_end = seed_next->read_start - 1;

      genome_start = seed_prev->genome_end + 1;
      genome_end = seed_prev->genome_end + read_gap + 5;      
      genome_read_sequence_by_chr_index(reference_prev, 0, 
					cal_next->chromosome_id - 1,
					&genome_start, &genome_end,
					genome);

      genome_start = seed_next->genome_start - read_gap - 5;
      genome_end = seed_next->genome_start - 1;
      genome_read_sequence_by_chr_index(reference_next, 0, 
					cal_next->chromosome_id - 1,
					&genome_start, &genome_end,
					genome);
      
      cigar_code_t *cigar_code_next = cal_next->info;

      //printf("REF-PREV: %s\n", reference_prev);
      //printf("REF-NEXT: %s\n", reference_next);
      //printf("read_gap = %i\n", read_gap);

      int c = 0;
      int r_pos = strlen(reference_next) - 1;

      distance = 0;
      while (distance < max_distance && 
	     c < read_gap) { 
	if (query_map[read_start + c] != reference_prev[c]) {	  
	  distance++;
	}
	c++;
      }
      prev_M = c;
      if (prev_M > 0) {
	cigar_code_append_new_op(prev_M, 'M', cigar_code_n);
      }
      //printf("Distance %i and Matches PREV %i\n", distance, prev_M);

      distance = 0;
      while (distance < max_distance && 
	     c < read_gap) {
	if (query_map[read_end - c] != reference_next[r_pos - c]) {	  
	  distance++;
	}
	c++;
      }
      next_M = c - prev_M;
      
      int num_I = read_end - (read_start + prev_M + next_M) + 1;

      genome_start = seed_prev->genome_end + prev_M;
      genome_end = seed_next->genome_start - next_M;

      if (num_I > 0) {
	cigar_code_append_new_op(num_I, 'I', cigar_code_n);
      }

      cigar_code_append_new_op(genome_end - genome_start - 1 - num_I, 'N', cigar_code_n);
            
      if (next_M > 0) {
	cigar_code_append_new_op(next_M, 'M', cigar_code_n);
      }
      //printf("Distance %i and Matches NEXT %i\n", distance, next_M);

      for (int i = 0; i < array_list_size(cigar_code_next->ops); i++) {
	cigar_op_t *op = array_list_get(i, cigar_code_next->ops);
	//printf("NEXT:: %i%c\n", op->number, op->name);
	cigar_code_append_new_op(op->number, op->name, cigar_code_n);
      }
            
      cal_prev = cal_next;
    }

  } else {
    return 1;
  }

  if (meta_alignment->cigar_right != NULL) {
    cigar_code_t *cigar_code_aux = meta_alignment->cigar_right;
    for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
      cigar_op_t *op = array_list_get(i, cigar_code_aux->ops);
      //printf("R:: %i%c\n", op->number, op->name);
      cigar_code_append_new_op(op->number, op->name, cigar_code_n);
    } 
    cigar_code_n->distance += cigar_code_aux->distance;
  } 
  
  meta_alignment->cigar_code = cigar_code_n;

  //printf("--->: %s\n", new_cigar_code_string(meta_alignment->cigar_code));
  return 0;
}

void meta_alignment_complete_free(meta_alignment_t *meta_alignment) {

  cigar_code_t *cigar_code = meta_alignment->cigar_code;
  seed_region_t *s;

  if (cigar_code != NULL) {
    for (int o = 0; o < array_list_size(cigar_code->ops); o++) {
      cigar_op_t *op = array_list_get(o, cigar_code->ops);
      if (op != NULL) cigar_op_free(op);
    }
  }
  
  for (int c = 0; c < array_list_size(meta_alignment->cals_list); c++) {
    cal_t *cal = array_list_get(c, meta_alignment->cals_list);	  
    //linked_list_item_t *item_list = cal->sr_list->first, *item_prev = NULL;
    while (s = linked_list_remove_first(cal->sr_list)) {
      if (s->info != NULL) {
	cigar_code = s->info;
	array_list_clear(cigar_code->ops, (void *)cigar_op_free);
	cigar_code_free(cigar_code);
	s->info = NULL;
      }
      seed_region_free(s);
    }
    linked_list_free(cal->sr_list, NULL);
    cal->sr_list = NULL;
    if (cal->info != NULL) { 
      cigar_code = cal->info;
      array_list_clear(cigar_code->ops, (void *)cigar_op_free);
      cigar_code_free(cigar_code); 
    }
    cal_free(cal);
  }
  array_list_free(meta_alignment->cals_list, NULL);
	  
  if (meta_alignment->cigar_left != NULL) {
    cigar_code = meta_alignment->cigar_left;
    array_list_clear(cigar_code->ops, (void *)cigar_op_free);
    cigar_code_free(cigar_code);
  }

  if (meta_alignment->cigar_right != NULL) {
    cigar_code = meta_alignment->cigar_right;
    array_list_clear(cigar_code->ops, (void *)cigar_op_free);
    cigar_code_free(cigar_code);
  }

  for (int c = 0; c < meta_alignment->num_cigars; c++) {
    cigar_code = meta_alignment->middle_cigars[c];
    if (cigar_code != NULL) {
      cigar_code = meta_alignment->middle_cigars[c];
      array_list_clear(cigar_code->ops, (void *)cigar_op_free);
      cigar_code_free(cigar_code); 
    }
  }

  meta_alignment_free(meta_alignment);

}

//NEW FUNCTION!
int apply_sw_rna(sw_server_input_t* input_p, batch_t *batch) {

  LOG_DEBUG("========= SPLICE JUNCTION SEARCH =========\n");

  int min_intron_size = input_p->min_intron_size;
  size_t max_intron_size = input_p->max_intron_size;
  
  //printf("%i - %i\n", min_intron_size, max_intron_size);
  
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  sw_optarg_t *sw_optarg = &input_p->sw_optarg;
  genome_t *genome = input_p->genome_p;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  size_t num_targets = mapping_batch->num_targets;
  metaexons_t *metaexons = input_p->metaexons;
  cal_optarg_t *cal_optarg = input_p->cal_optarg_p;
  avls_list_t *avls_list = input_p->avls_list;
  bwt_optarg_t *bwt_optarg = input_p->bwt_optarg_p;
  bwt_index_t *bwt_index = input_p->bwt_index_p;

  //if (!bwt_index->dirname) { exit(-1); }
  linked_list_t *buffer = input_p->buffer;
  linked_list_t *buffer_hc = input_p->buffer_hc;
  FILE *f_sa = input_p->f_sa;
  FILE *f_hc = input_p->f_hc;

  array_list_t *cals_list, *fusion_cals, *fusion_cals_aux;
  cal_t *cal, *cal_prev, *cal_next, *first_cal, *last_cal;
  fastq_read_t *fq_read;
  linked_list_iterator_t itr;  
  seed_region_t *s, *s_prev, *s_next;
  cigar_code_t *cigar_code, *cigar_code_prev, *cigar_code_aux;
  cigar_code_t *alig_cigar_code;
  cigar_op_t *cigar_op_start, *cigar_op_end, *cigar_op, *cigar_op_prev, *cigar_op_aux;
  int *new_targets = (int *)calloc(mapping_batch->num_allocated_targets, sizeof(int));
  array_list_t *merge_cals;
  linked_list_t *linked_list;
  seed_region_t *seed_region;

  //float *cals_score = (float *)calloc(100, sizeof(float));
  float score;
  char reference[2048];
  //char reference_prev[2048];
  //char reference_next[2048];
  //char reference_aux[2048];
  char query[2048];
  //char query_revcomp[2048];
  alignment_t *alignment;
  char q[2048];
  char r[2048];

  //char **rev_comp = (char **)calloc(num_reads, sizeof(char *));

  char *sequence;
  char *query_ref;
  char *quality_map, *query_map;
  //float scores_ranking[mapping_batch->num_allocated_targets][50];
  float *cals_score;
  //char cigar_str[1024];

  cigar_op_t *first_op;
  char *match_seq, *match_qual, *optional_fields, *p;
  int match_start, match_len, optional_fields_length, AS;
  float norm_score;

  int delete_not_cigar;
  int padding_left, padding_right, len_query;  
  int num_sw = 0;
  int num_sw_sp = 0;
  int num_sw_ext = 0;
  size_t num_new_targets = 0;
  int coverage;
  int read_length;
  int cal_pos = 0;
  int number_of_best;

  //int num_extrem_ok = 0;
  //int num_sp_ok = 0;

  //Change to report most alignments
  int seed_err_size = 20;
  int n_alignments = 1;
  int target;
  int c, exact_nt;
  size_t genome_start, genome_end, genome_start2, genome_end2, read_start, read_end;
  int flank = 20;
  int n_delete;
  int sw_pos;
  int seeds_nt; 
  int e1;
  int e2;
  int ref_pos;
  int nt_2;
  
  int start_seeding, end_seeding;
  int lim_start, lim_end;

  metaexon_t *first_metaexon, *last_metaexon;
  int l_found, r_found; 
  size_t num_cals;
  int t;
  register int i, j;

  int meta_type;
  meta_alignment_t *meta_alignment;
  sw_item_t *sw_item;  
  sw_depth_t sw_depth;
  sw_depth.depth = 0;
  sw_multi_output_t *output = sw_multi_output_new(MAX_DEPTH);

  array_list_t **meta_alignments_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  array_list_t *seeds_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *final_positions = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED); 
  int make_seeds = 0;

  array_list_t *cals_targets = array_list_new(500, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  int *post_process_reads = (int *)calloc(num_reads, sizeof(int));
  int *data_type  = (int *)calloc(num_reads, sizeof(int));
  //array_list_t **data_file_reads = (array_list_t **)calloc(num_reads, sizeof(array_list_t *));

  float **scores_ranking = (float **)calloc(num_reads, sizeof(float *));//[num_reads][200];
  int read_nt;

  int seed_size = cal_optarg->seed_size;

  int pair_mode = input_p->pair_mode;  
  // array_list flag: 0 -> Not  BWT Anchors found (NOT_ANCHORS)        *
  //                  1 -> One  BWT Anchors found (SINGLE_ANCHORS)     *
  //                  2 -> Pair BWT Anchors found (DOUBLE_ANCHORS)     *
  //                  3 -> Alignments found       (ALIGNMENTS_FOUND)
  //                  4 -> Alignments exceded     (ALIGNMENTS_EXCEEDED)
  int flag;

  //pthread_mutex_lock(&mutex_sp);
  //TOTAL_READS_PROCESS += num_reads;
  //pthread_mutex_unlock(&mutex_sp);



  for (i = 0; i < num_reads; i++) {        
    cals_list = mapping_batch->mapping_lists[i];
    flag = array_list_get_flag(cals_list);
    fq_read = array_list_get(i, mapping_batch->fq_batch);

    //printf("MODE: %i\n", flag);

    num_cals = array_list_size(cals_list);
    meta_alignments_list[i] = array_list_new(num_cals,
				     1.25f,
				     COLLECTION_MODE_ASYNCHRONIZED);

    scores_ranking[i] = (float *)calloc(num_cals + 10, sizeof(float));//[num_reads][200];

    //printf("%s : flag %i : %i\n", fq_read->id, flag, array_list_size(cals_list));
    //if (flag == ALIGNMENTS_FOUND || flag == ALIGNMENTS_EXCEEDED) {
    //data_type[i] = ALIGNMENT_TYPE;
    if (flag == DOUBLE_ANCHORS) {
      //printf("\tWK_1ph: -- DOUBLE ANCHOR PROOCESS --\n");
      //printf("<<<<@ %s\n", fq_read->id);
      //1st- Merge Double anchors CALs. TODO: make only merge, not generate score	
      //array_list_clear(mapping_batch->mapping_lists[i], cal_free);continue;
      make_seeds = 0;
      int found[array_list_size(cals_list)];
      int founds = 0;
      for (int p = 0; p < array_list_size(cals_list); p += 2) {	
	first_cal = array_list_get(p, cals_list);
	//cal_print(first_cal);
	last_cal  = array_list_get(p + 1, cals_list);
	//cal_print(last_cal);
	if (first_cal->end + min_intron_size/*fq_read->length*/ >= last_cal->start) {
	  found[p] = 1;
	  found[p + 1] = 1;
	  founds++;
	} else {
	  found[p] = 0;
	  found[p + 1] = 0;	  
	}
      }

      if (founds) {
	//array_list_t *cals_aux = array_list_new(array_list_size(cals_list), 1.5f, COLLECTION_MODE_ASYNCHRONIZED);
	//printf("REPORT FILL GAPS\n");
	for (int p = 0; p < array_list_size(cals_list); p += 2) {
	  first_cal = array_list_get(p, cals_list);
	  if (first_cal->strand == 1) {
	    //if (!rev_comp[i]) {
	    //rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	    //strcpy(rev_comp[i], fq_read->sequence);
	    //seq_reverse_complementary(rev_comp[i], fq_read->length);
	    //}
	    query_map = fq_read->revcomp;
	  } else {
	    query_map = fq_read->sequence;
	  }
	  last_cal = array_list_get(p + 1, cals_list);
	  if (found[p] == 1) {
	    s = linked_list_get_first(last_cal->sr_list);
	    seed_region = seed_region_new(s->read_start, s->read_end, 
					  s->genome_start, s->genome_end, 1, 0, 0);
	    linked_list_insert_last(seed_region, first_cal->sr_list);
	    first_cal->end = last_cal->end;
	    cal_free(last_cal);
	    meta_alignment = meta_alignment_cal_new(first_cal);
	    meta_alignment_fill_gaps(META_ALIGNMENT_NONE,
				     meta_alignment, query_map, genome,
				     sw_optarg, output, metaexons, 
				     &sw_depth, avls_list, min_intron_size);	    
	    array_list_insert(meta_alignment, meta_alignments_list[i]);
	  } else {
	    cal_free(first_cal);
	    cal_free(last_cal);
	  }
	}
	continue;
      } else {
	first_cal = array_list_get(0, cals_list);
	last_cal  = array_list_get(1, cals_list);
	s_prev = linked_list_get_last(first_cal->sr_list);    
	s_next = linked_list_get_first(last_cal->sr_list);	
	int gap = s_next->read_start - s_prev->read_end;
	if (gap > 20) {
	  make_seeds = 1;
	}
      }

      int found_sp;
      int first_cal_len, last_cal_len;
      //array_list_t *delete_targets = array_list_new(array_list_size(cals_list), 
      //					    1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    
      for (int p = 0; p < array_list_size(cals_list); p += 2) {
	first_cal = array_list_get(p, cals_list);
	first_cal_len = first_cal->end - first_cal->start;
	last_cal  = array_list_get(p + 1, cals_list);
	last_cal_len = last_cal->end - last_cal->start;
	if (last_cal->strand == 1) {
	  //if (!rev_comp[i]) {
	  //rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	  //strcpy(rev_comp[i], fq_read->sequence);
	  //seq_reverse_complementary(rev_comp[i], fq_read->length);
	  //}
	  query_map = fq_read->revcomp;
	  //query_map = rev_comp[i];
	} else {
	  query_map = fq_read->sequence;
	}
	/*
	printf("--CALs-- Iter %i --\n", p);
	cal_print(first_cal);
	cal_print(last_cal);
	printf("--CALs-- Iter %i --\n", p);
	*/
	s_prev = linked_list_get_last(first_cal->sr_list);    
	s_next = linked_list_get_first(last_cal->sr_list);
	
	cigar_code = search_double_anchors_cal(query_map,
					       first_cal, last_cal,
					       metaexons, genome,
					       fq_read, &meta_type, 
					       avls_list);
	
	//cigar_code = NULL;
	if (cigar_code != NULL) { 
	  meta_alignment = meta_alignment_new();
	  meta_alignment_insert_cal(first_cal, meta_alignment);
	  if (meta_type == META_ALIGNMENT_LEFT) {
	    meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
	    cal_free(last_cal);
	  } else {
	    meta_alignment_insert_cal(last_cal, meta_alignment);
	    meta_alignment_insert_cigar(cigar_code, CIGAR_SIMPLE_MIDDLE, 0, meta_alignment);
	  }
	  meta_alignment_fill_gaps(meta_type,
				   meta_alignment, query_map, genome,
				   sw_optarg, output, metaexons, 
				   &sw_depth, avls_list, min_intron_size);	  
	  array_list_insert(meta_alignment, meta_alignments_list[i]);
	} else {
	  //====(((((S M I T H - W A T E R M A N    P H A S E)))))====//	
	  //TODO: Extend before caling
	  cal_t *father_cal;
	  int gap = s_next->read_start - s_prev->read_end;
	  array_list_t *init_list = array_list_new(40, 1.25f, 
						   COLLECTION_MODE_ASYNCHRONIZED);
	  array_list_t *aux_list = array_list_new(40, 1.25f, 
						  COLLECTION_MODE_ASYNCHRONIZED);
	  region_t *region_prev, *region_next;
	  if (make_seeds) {
	    int read_start, read_end;
	    if (first_cal->strand == 1) {
	      read_start = fq_read->length - s_next->read_start;
	      read_end = (fq_read->length - (s_prev->read_end + 1)) - 1;
	    } else {
	      read_start = s_prev->read_end + 1;
	      read_end = s_next->read_start - 1;
	    }
	    //printf("SEEDING BETWEEN +(%i - %i), -(%i - %i)\n", s_prev->read_end + 1, s_next->read_start - 1, read_start, read_end);
	    bwt_map_exact_seeds_by_region(read_start, read_end,
					  fq_read->sequence, 16, 16,
					  bwt_optarg, bwt_index,
					  seeds_list);
	    make_seeds = 0;
	  }

	  //array_list_t *mapping_list = array_list_new(array_list_size(seeds_list) + 2, 1.25f, 
	  //COLLECTION_MODE_ASYNCHRONIZED);
	  if (array_list_size(seeds_list) > 0) {
	    generate_cals_between_anchors(CALING_DOUBLE_ANCHORS,
					  genome->num_chromosomes,
					  first_cal, 
					  last_cal, 
					  fq_read, 
					  seeds_list, 
					  aux_list, 
					  cal_optarg);
	    
	    if (array_list_size(aux_list) > 5) {
	      for (int zz = array_list_size(aux_list) - 2; zz > 0; zz--) {
		cal_t *cal_aux = array_list_remove_at(zz, aux_list);
		cal_free(cal_aux);
	      }
	    }
	  } else {
	    array_list_insert(first_cal, aux_list);
	    array_list_insert(last_cal, aux_list);
	  }
	  /*
	  for (int zz = 0; zz < array_list_size(aux_list); zz++) {
	    cal_print(array_list_get(zz, aux_list));
	  }
	  */
	  if (array_list_size(aux_list) > 0) {
	    first_cal = array_list_get(0, aux_list);
	    //first_cal->fill = 1;
	    if (array_list_size(aux_list) == 1) {
	      meta_alignment = meta_alignment_cal_new(first_cal);
	      //meta_alignment_set_status(META_CLOSE, meta_alignment);
	      meta_type = META_ALIGNMENT_NONE;
	    } else {
	      meta_alignment = meta_alignment_cals_new(aux_list);		
	      meta_type = META_ALIGNMENT_MIDDLE;
	      for (int c = 1; c < array_list_size(aux_list); c++) { 
		last_cal = array_list_get(c, aux_list);

		avl_node_t *node_avl_start, *node_avl_end;
		seed_region_t *s_prev_aux = linked_list_get_last(first_cal->sr_list);    
		seed_region_t *s_next_aux = linked_list_get_first(last_cal->sr_list);
		int distance_aux;

		size_t sp_start, sp_end;
		int sp_type;
		int nt;

		nt = search_simple_splice_junction(s_prev_aux, s_next_aux,
						   first_cal->chromosome_id, 
						   first_cal->strand,
						   query_map, genome, 
						   &sp_start, &sp_end,
						   &sp_type,
						   &distance_aux);

		//printf("nt = %i, %lu < %lu, %lu < %lu\n", nt, s_prev_aux->genome_start, node_avl_start->position - 1, 
		//     node_avl_end->position + 1, s_next_aux->genome_end);
		//nt = 0;
		if (nt &&
		    s_prev_aux->genome_start < sp_start &&
		    sp_end < s_next_aux->genome_end) {
		  //s_prev_aux->genome_start < node_avl_start->position - 1 &&
		  //node_avl_end->position + 1 < s_next_aux->genome_end) {
		  
		  cigar_code = cigar_code_new();
		  cigar_code->distance = distance_aux;
		  cigar_code_append_new_op(nt, 'N', cigar_code);
		  int err_l, err_r;

		  //assert(s_prev_aux->genome_start < sp_start);
		  //assert(node_avl_end->position + 1 < s_next_aux->genome_end);
		  int sp_strand = (sp_type == CT_AC_SPLICE ? 1 : 0);

		  allocate_start_node(first_cal->chromosome_id - 1,
				      sp_strand,
				      sp_start,
				      sp_end,
				      sp_start,
				      sp_end,
				      FROM_READ,
				      sp_type,
				      NULL, 
				      &node_avl_start,
				      &node_avl_end, 
				      avls_list);
		  
		  err_l = metaexon_insert(1, first_cal->chromosome_id - 1,
					  s_prev_aux->genome_start,
					  node_avl_start->position - 1,
					  min_intron_size,
					  METAEXON_RIGHT_END, node_avl_start,
					  metaexons);
		  
		  err_r = metaexon_insert(1, last_cal->chromosome_id - 1,
					  node_avl_end->position + 1,
					  s_next_aux->genome_end, 
					  min_intron_size,
					  METAEXON_LEFT_END, node_avl_end,
					  metaexons);
		  
		  meta_alignment_insert_cigar(cigar_code, CIGAR_SIMPLE_MIDDLE, c - 1, meta_alignment);

		} else { 			  
		  //pthread_mutex_lock(&mutex_sp);
		  //printf(":@ Not Found:\n");
		  //pthread_mutex_unlock(&mutex_sp);
		  
		  //printf("2.PROCESS.2 $ %s\n", fq_read->id);
		  
		  info_sp_t *info_sp = sw_reference_splice_junction(first_cal, last_cal,
								    query_map, genome,
								    q, r);
		  //New sw item. Storing data...
		  sw_item = sw_item_new(SP_SW, i, 0, c - 1,
					first_cal, last_cal, 
					meta_alignment, NULL, 
					NULL, info_sp);	      
		  //Insert item... and process if depth is full
		  sw_depth_insert(q, r, sw_item,
				  sw_optarg, output,
				  avls_list, metaexons, 
				  &sw_depth, genome, min_intron_size);	
		  //printf("\t-->(%i) SW\n", sw_depth.depth);
		  //meta_alignment_insert_cigar(NULL, CIGAR_SIMPLE_MIDDLE, c - 1, meta_alignment);
		}

		first_cal = last_cal;
	      } //End loop aux list 
	    }
	    meta_alignment_fill_gaps(meta_type,
				     meta_alignment, query_map, genome,
				     sw_optarg, output, metaexons, 
				     &sw_depth, avls_list, min_intron_size);
	    array_list_insert(meta_alignment, meta_alignments_list[i]);
	  }
	  array_list_free(aux_list, NULL);
	  array_list_free(init_list, NULL);
	} //Not found splice in metaexon      
	//====(((((S M I T H - W A T E R M A N    P H A S E    E N D)))))====//
      } //End loop double anchors (cals_list)      

      if (array_list_size(meta_alignments_list[i]) == 0) {
	array_list_clear(mapping_batch->mapping_lists[i], (void *)NULL);
      }

    } else if (flag == SINGLE_ANCHORS) {
      //array_list_clear(cals_list, (void*)cal_free);
      //continue;
      //printf("\tWK_1ph: -- SINGLE ANCHOR PROOCESS --\n");
      int map;      
      int read_nt;
      int seeds_process = 0;
      cigar_code_t *cigar_code_res;
      array_list_t *delete_targets = array_list_new(array_list_size(cals_list), 1.25f,
						    COLLECTION_MODE_ASYNCHRONIZED);

      for (size_t p = 0; p < array_list_size(cals_list); p++) {
	cal = array_list_get(p, cals_list);	
	//cal_print(cal);
	//map = 0;
	cigar_code_res = NULL;
	meta_alignment = NULL;
	if (cal->strand == 1) {
	  //if (!rev_comp[i]) {
	  //rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	  //strcpy(rev_comp[i], fq_read->sequence);
	  //seq_reverse_complementary(rev_comp[i], fq_read->length);
	  //}
	  query_map = fq_read->revcomp;
	  //query_map = rev_comp[i];
	} else {
	  query_map = fq_read->sequence;
	}
	
	s_prev = linked_list_get_last(cal->sr_list);
	if (s_prev->read_start == 0) {
	  read_nt = fq_read->length - (s_prev->read_end + 1);
	} else {
	  read_nt = s_prev->read_start;
	}
		
	if (s_prev->read_start == 0) {
	  cigar_code_res = fill_extrem_gap(query_map, 
					   cal,
					   FILL_GAP_RIGHT,
					   genome,
					   metaexons, 
					   avls_list);
	} else {
	  cigar_code_res = fill_extrem_gap(query_map, 
					   cal,
					   FILL_GAP_LEFT,
					   genome,
					   metaexons, 
					   avls_list);	
	}

	if (cigar_code_res) {
	  meta_alignment = meta_alignment_new();
	  meta_alignment_insert_cal(cal, meta_alignment);
	  
	  if (s_prev->read_start == 0) {
	    meta_alignment->cigar_right = cigar_code_res;
	  } else {
	    meta_alignment->cigar_left = cigar_code_res;
	  }
	  
	  meta_alignment_fill_gaps(META_ALIGNMENT_RIGHT,
				   meta_alignment, query_map, genome,
				   sw_optarg, output, metaexons, 
				   &sw_depth, avls_list, min_intron_size);
	  
	  array_list_insert((void *)meta_alignment, meta_alignments_list[i]);
	}

	if (!meta_alignment) { array_list_insert((void *)p, delete_targets); } 

      }
	
      //---------------------------------//
      
      if (array_list_size(meta_alignments_list[i]) <= 0) { 
	post_process_reads[i] = 1;
	data_type[i] = CAL_TYPE;
	//pthread_mutex_lock(&mutex_sp);
	//insert_file_item(fq_read, cals_list, f_sa);
	//pthread_mutex_unlock(&mutex_sp);

	//array_list_clear(cals_targets, (void*)NULL);	
	//pthread_mutex_lock(&mutex_sp);
	//extern size_t w2_r;
	//w2_r++;
	//file_write_items(fq_read, mapping_batch->mapping_lists[i],
	//		 CAL_TYPE, f_sa, f_hc, 0);
	//pthread_mutex_unlock(&mutex_sp);
	
	//No clear ??
	//array_list_clear(mapping_batch->mapping_lists[i], (void*)cal_free);
	
	//pthread_mutex_lock(&mutex_sp);
	//TOTAL_READS_SA++;
	//pthread_mutex_unlock(&mutex_sp);

      } else {
	for (t = 0; t < array_list_size(delete_targets); t++) {
	  size_t target = (size_t)array_list_get(t, delete_targets);
	  cal_t *cal = (cal_t *)array_list_get(target, cals_list);
	  cal_free(cal);
	  array_list_set(target, NULL, cals_list);
	}
      }

      array_list_free(delete_targets, NULL);

    } else if (flag == NOT_ANCHORS) {
      //printf( "\tWK_1ph: -- NOT ANCHORS (CALS PROCESS) --\n");
      //printf("NUM CALs %i\n", array_list_size(cals_list));
      //array_list_clear(cals_list, cal_free);
      //continue;
      if (array_list_size(cals_list) <= 0) {
	continue;
      }

      merge_and_filter_cals(cals_targets, cals_list, fq_read, 
			    bwt_optarg, bwt_index, genome, 
			    scores_ranking[i]);

      int MAX_PROCESS = 2;
      float best_score = scores_ranking[i][0];

      for (int j = 1; j < array_list_size(cals_targets); j++) { 
	if (best_score == scores_ranking[i][j]) {
	  MAX_PROCESS++;
	}
      }

      if (MAX_PROCESS > 5) { MAX_PROCESS = 5; }

      //printf("NOW NUM CALs %i\n", cals_targets->size);
      if (array_list_size(cals_targets) > MAX_PROCESS) {      
	//size_t seed_size = 16;
	//float best_score = scores_ranking[i][0];
	//int limit;	
	//int reorder = 0;

	/*if (best_score < 40.0) {
	  for (int j = 0; j < array_list_size(cals_targets); j++) {	   
	    //array_list_t *fusion_list = array_list_get(j, cals_targets);
	    //cal = array_list_get(0, fusion_list);
	    //seed_region_t *s_prev = linked_list_get_first(cal->sr_list);    
	    //seed_region_t *s_next = linked_list_get_last(cal->sr_list);	
	    //int last_nt = fq_read->length % seed_size;
	    //int cal_size = s_next->read_end - s_prev->read_start;
	    /*if (s_prev->read_start == (seed_size / 2) || 
		s_prev->read_end == ((seed_size / 2) + seed_size) || 
		s_prev->read_start == last_nt) {
	      scores_ranking[i][j] = 0;
	      reorder = 1;
	    } else if (last_nt != 0 &&
		       s_next->read_end == (fq_read->length - last_nt - 1) || 
		       s_next->read_start == (read_length - last_nt - seed_size)) {
	      scores_ranking[i][j] = 0;
	      reorder = 1;	    
	      }
	    scores_ranking[i][j] >= seed_size*2)
	  }
	}

	if (reorder) {
	  float score;
	  for (int z = 0; z < array_list_size(cals_targets) - 1; z++) {
	    for (j = i + 1; j < array_list_size(cals_targets); j++) {
	      if (scores_ranking[i][j] > scores_ranking[i][z]) {
		array_list_swap(z, j, cals_targets);
		score = scores_ranking[i][j];
		scores_ranking[i][j] = scores_ranking[i][z];
		scores_ranking[i][z] = score;
	      }
	    }
	  }
	  }*/
      
	for (int j = array_list_size(cals_targets) - 1; j >= MAX_PROCESS; j--) { 
	  array_list_t *fusion_list = array_list_remove_at(j, cals_targets);
	  array_list_free(fusion_list, (void *)cal_free);
	}	

      }     

      //printf("::===:: BEST SCORE %f WITH %i CALs ::===:: %s\n", scores_ranking[i][0], 
      //   array_list_size(cals_targets), fq_read->id);      

      int meta_type;


      //order_cals(cals_list);
      int limit = array_list_size(cals_targets);
      //MAX_PROCESS > array_list_size(cals_targets) ?
      //array_list_size(cals_targets) : MAX_PROCESS;	
      //printf("::----> Process LIM %i\n", limit);
      for (int j = 0; j < limit; j++) {
	meta_alignment = NULL;
	meta_type = -1;
	array_list_t *fusion_list = array_list_get(j, cals_targets);
	/*
	printf("==== FUSION CALS PROCESS ====\n");
	for (int t = 0; t < array_list_size(fusion_list); t++) {
	  cal_t *cal_aux = array_list_get(t, fusion_list);
	  cal_print(cal_aux);
	}
	printf("==== ------------------- ====\n");
	*/
	first_cal = array_list_get(0, fusion_list); 
	if (first_cal->strand == 1) {
	  //if (!rev_comp[i]) {
	  //rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	  //strcpy(rev_comp[i], fq_read->sequence);
	  //seq_reverse_complementary(rev_comp[i], fq_read->length);
	  //}
	  //query_map = rev_comp[i];
	  query_map = fq_read->revcomp;
	} else {
	  query_map = fq_read->sequence;
	}
	meta_alignment = NULL;
	if (array_list_size(fusion_list) > 1) { 
	  //printf("FUSION CALS REPORT\n");
	  cal_t *first_cal = array_list_get(0, fusion_list);
	  cal_t *last_cal  = array_list_get(array_list_size(fusion_list) - 1, fusion_list);
	  
	  //TODO: IF WE HAVE MORE THAN TWO CALS SEARCH SINGLE ANCHOR
	  cigar_code = search_double_anchors_cal(query_map,
						 first_cal, last_cal,
						 metaexons, genome,
						 fq_read, &meta_type, avls_list);

	  //cigar_code = NULL;
	  if (cigar_code != NULL) {
	    meta_alignment = meta_alignment_new();
	    meta_alignment_insert_cal(first_cal, meta_alignment);
	    if (meta_type == META_ALIGNMENT_LEFT) {
	      meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
	      for (int c = 1; c < array_list_size(fusion_list); c++) {
		cal = array_list_get(c, fusion_list);
		cal_free(cal);
		array_list_set(c, NULL, fusion_list);
	      }
	      fusion_list->size = 1;
	    } else {
	      for (int c = 1; c < array_list_size(fusion_list); c++) {
		cal = array_list_get(c, fusion_list);
		meta_alignment_insert_cal(cal, meta_alignment);
	      }
	      meta_alignment_insert_cigar(cigar_code, CIGAR_SIMPLE_MIDDLE, 0, meta_alignment);
	    }
	    //array_list_free(fusion_list, NULL);
	    meta_alignment_fill_gaps(meta_type,
				     meta_alignment, query_map, genome,
				     sw_optarg, output, metaexons, 
				     &sw_depth, avls_list, min_intron_size);	
	  } else {
	    //printf("NOT FOUND\n");
	    meta_alignment = meta_alignment_cals_new(fusion_list);  
	    for (int c = 1; c < array_list_size(fusion_list); c++) { 
	      last_cal = array_list_get(c, fusion_list);
	      avl_node_t *node_avl_start, *node_avl_end;
	      seed_region_t *s_prev = linked_list_get_last(first_cal->sr_list);    
	      seed_region_t *s_next = linked_list_get_first(last_cal->sr_list);	
	      int distance_aux;
	      size_t sp_start, sp_end;
	      int sp_type;
	      int nt = search_simple_splice_junction(s_prev, s_next, 
						     last_cal->chromosome_id, 
						     last_cal->strand,
						     query_map, genome, 
						     &sp_start, &sp_end,
						     &sp_type,
						     &distance_aux);

	      if (nt) {
		//printf("Found Splice simple\n");
		cigar_code = cigar_code_new();
		cigar_code->distance = distance_aux;
		cigar_code_append_new_op(nt, 'N', cigar_code);
		
		seed_region_t *s_prev_aux = linked_list_get_last(first_cal->sr_list);    
		seed_region_t *s_next_aux = linked_list_get_first(last_cal->sr_list);
		
		int err_l, err_r;

		int sp_strand = (sp_type == CT_AC_SPLICE ? 1 : 0);

		allocate_start_node(first_cal->chromosome_id - 1,
				    sp_strand,
				    sp_start,
				    sp_end,
				    sp_start,
				    sp_end,
				    FROM_READ,
				    sp_type,
				    NULL, 
				    &node_avl_start,
				    &node_avl_end, 
				    avls_list);		  
		
		assert(s_prev_aux->genome_start < node_avl_start->position - 1);
		assert(node_avl_end->position + 1 < s_next_aux->genome_end);
		
		err_l = metaexon_insert(first_cal->strand, first_cal->chromosome_id - 1,
					s_prev_aux->genome_start, node_avl_start->position - 1, min_intron_size,
					METAEXON_RIGHT_END, node_avl_start,
					metaexons);
		
		err_r = metaexon_insert(last_cal->strand, last_cal->chromosome_id - 1,
					node_avl_end->position + 1, s_next_aux->genome_end, min_intron_size,
					METAEXON_LEFT_END, node_avl_end,
					metaexons);
		
		meta_alignment_insert_cigar(cigar_code, CIGAR_SIMPLE_MIDDLE, c - 1, meta_alignment);
		
	      } else { 
		//pthread_mutex_lock(&mutex_sp);
		//printf(":@ Not Found:\n");
		//pthread_mutex_unlock(&mutex_sp);
		
		
		info_sp_t *info_sp = sw_reference_splice_junction(first_cal, last_cal,
								  query_map, genome,
								  q, r);
		//New sw item. Storing data...
		sw_item = sw_item_new(SP_SW, i, 0, c - 1,
				      first_cal, last_cal, 
				      meta_alignment, NULL, 
				      NULL, info_sp);	      
		
		//Insert item... and process if depth is full
		sw_depth_insert(q, r, sw_item,
				sw_optarg, output,
				avls_list, metaexons,
				&sw_depth, genome, min_intron_size);
		
		//meta_alignment_insert_cigar(NULL, CIGAR_SIMPLE_MIDDLE, c - 1, meta_alignment);
		//printf("\t-->(%i) SW\n", sw_depth.depth);
	      }
	      first_cal = last_cal;
	    }
	    //array_list_free(fusion_cals, NULL);
	    meta_alignment_fill_gaps(META_ALIGNMENT_MIDDLE,
				     meta_alignment, query_map, genome,
				     sw_optarg, output, metaexons, 
				     &sw_depth, avls_list, min_intron_size);	 
	  }
	  
	  if (meta_alignment) {
	    array_list_insert(meta_alignment, meta_alignments_list[i]);
	  }
	} else {
	  //printf(":::: SINGLE MAP (%i)::::\n", array_list_size(fusion_list));
	  cal = array_list_get(0, fusion_list);
	  //cal_print(cal);
	  
	  seed_region_t *s_prev = linked_list_get_first(cal->sr_list);
	  seed_region_t *s_next = linked_list_get_last(cal->sr_list);
	  int gap_start = s_prev->read_start;
	  int gap_end = s_next->read_end - fq_read->length - 1;	      
	  int process_left = 0, process_right = 0;

	  meta_alignment = meta_alignment_new();
	  meta_alignment_insert_cal(cal, meta_alignment);
	  meta_alignment_fill_gaps(META_ALIGNMENT_MIDDLE,
				   meta_alignment, query_map, genome,
				   sw_optarg, output, metaexons, 
				   &sw_depth, avls_list, min_intron_size);	 	      
	  array_list_insert(meta_alignment, meta_alignments_list[i]);		
	}

	cal_t *first_cal = array_list_get(0, fusion_list);
	cal_t *last_cal  = array_list_get(array_list_size(fusion_list) - 1, fusion_list);	    
	
	assert(first_cal != NULL);
	assert(last_cal != NULL);
	seed_region_t *s_prev = linked_list_get_first(first_cal->sr_list);
	seed_region_t *s_next = linked_list_get_last(last_cal->sr_list);
	cigar_code_t *cc_left, *cc_right;
	
	if (s_prev->read_start != 0 && 
	    meta_alignment->cigar_left == NULL) {
	  cc_left = fill_extrem_gap(query_map, 
				    first_cal,
				    FILL_GAP_LEFT,
				    genome,
				    metaexons, 
				    avls_list); 
	  meta_alignment_insert_cigar(cc_left, CIGAR_ANCHOR_RIGHT, 0, meta_alignment);
	}
	
	if (meta_alignment->cigar_right == NULL &&
	    s_next->read_end != fq_read->length - 1) {
	  cc_right = fill_extrem_gap(query_map, 
				     last_cal,
				     FILL_GAP_RIGHT,
				     genome,
				     metaexons, 
				     avls_list); 
	  meta_alignment_insert_cigar(cc_right, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
	}

      }
    
      
      if (array_list_size(meta_alignments_list[i]) > 0) {
	for (int j = 0; j < array_list_size(cals_targets); j++) {
	  array_list_t *fusion_list = array_list_get(j, cals_targets);
	  array_list_free(fusion_list, NULL);
	}
      }
    
      array_list_clear(mapping_batch->mapping_lists[i], (void*)NULL);
      array_list_clear(cals_targets, (void*)NULL);
      
    }
    
    array_list_clear(seeds_list, (void*)region_bwt_free);
    
  } //End loop reads
  
  sw_depth_process(sw_optarg, output, 
		   &sw_depth, avls_list, metaexons, SW_FINAL, 
		   genome, min_intron_size);


  //Merge splice junctions & order by distance
  //printf("CLOSE META ALIGNMENTS\n");

  for (int i = 0; i < num_reads; i++) {
    fq_read = array_list_get(i, mapping_batch->fq_batch);
    //printf("<<<<CLOSE META (%i) %s:\n", array_list_size(meta_alignments_list[i]), fq_read->id);

    //delete else
    if (array_list_size(meta_alignments_list[i]) == 0) { 
      continue; 
    } 

    //printf("Num meta %i\n", array_list_size(meta_alignments_list[i]));
    for (int m = 0; m < array_list_size(meta_alignments_list[i]); m++) {
      meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);
      if (meta_alignment_get_status(meta_alignment) == META_OPEN) {
	meta_alignment_close(meta_alignment);     
      }
    }

    //TODO: Alert! Reads with middle small exons in middle no map!Add 'I' operation in the cigar
    int no_map = 1;
    for (int m = 0; m < array_list_size(meta_alignments_list[i]); m++) {
      meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);
      //printf("\t CIGAR: %s\n", new_cigar_code_string(meta_alignment->cigar_code));
      if (meta_alignment->score == fq_read->length) {
	no_map = 0;
	break;
      } 
    }
    
    meta_alignment_t *meta_alignment = array_list_get(0, meta_alignments_list[i]);
    if (no_map) {
      //printf("NO MAP\n");
      meta_alignment->flag = 0;
    } else {
      //printf("MAP\n");
      meta_alignment->flag = 1;
    }

  }
  
  //return RNA_POST_PAIR_STAGE;

  //printf("WK_1ph: =============== REPORT ALIGNMENTS =====================\n");
  //printf("REPORT META\n");
  size_t start_mapping;
  for (int i = 0; i < num_reads; i++) {
    fq_read = array_list_get(i, mapping_batch->fq_batch);

    //printf("(%i) WK_1ph-Meta_Report:  %s >>>\n", array_list_size(meta_alignments_list[i]), fq_read->id);
    
    if (array_list_size(meta_alignments_list[i]) > 0) {
      array_list_clear(mapping_batch->mapping_lists[i], (void*)NULL);
    } else {
      continue;
    }

    int map = 0;
    meta_alignment_t *meta_alignment = array_list_get(0, meta_alignments_list[i]);
    if (meta_alignment->flag == 1) {      
      data_type[i] = ALIGNMENT_TYPE;
      for (int m = 0; m < array_list_size(meta_alignments_list[i]); m++) {
	meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);
   	map = 1;
	char query[2048];
	char quality[2048];
	int map = 0;

	//if (meta_alignment_get_status(meta_alignment) == META_CLOSE) {	  
	optional_fields_length = 0;
	optional_fields = NULL;
	
	first_cal = meta_alignment_get_first_cal(meta_alignment);
	s_prev = linked_list_get_first(first_cal->sr_list);
	int h_left = s_prev->read_start;
	
	last_cal = meta_alignment_get_last_cal(meta_alignment);
	s_next = linked_list_get_last(last_cal->sr_list);
	int h_right = (fq_read->length - 1) - s_next->read_end;

	if (first_cal->strand == 1) {
	  //query_map = rev_comp[i];
	  query_map = fq_read->revcomp;
	} else {
	  query_map = fq_read->sequence;
	}

	cigar_code = meta_alignment->cigar_code;	
	assert(cigar_code != NULL);


	//printf("START_LEFT = %i, START_RIGHT = %i, start_mapping = %lu, ", h_left, h_right, first_cal->start);
	start_mapping = first_cal->start;
	int dsp = 0;

	if (h_left > 0 && 
	    meta_alignment->cigar_left == NULL) {
	  //printf("H_LEFT ->  %i -> seed %i\n", h_left, s_prev->read_start);
	  array_list_insert_at(0, cigar_op_new(h_left, 'H'), cigar_code->ops);
	} else {
	  h_left = 0;
	  if (meta_alignment->cigar_left != NULL ) {
	    cigar_code_t *cigar_code = meta_alignment->cigar_left;
	    for (int c = 0; c < cigar_code->ops->size; c++) {
	      cigar_op_t *op = array_list_get(c, cigar_code->ops);
	      if (op->name == 'M' ||
		  op->name == 'N' ||
		  op->name == 'D') {
		dsp += op->number;
	      }
	    }
	  }	 	  
	}
	
	start_mapping -= dsp;
	//printf(" new_start = %lu\n", start_mapping);

	if (h_right > 0 &&
	    meta_alignment->cigar_right  == NULL) {
	  //printf("H_RIGHT ->  %i -> seed %i\n", h_right,  s_next->read_end);
	  array_list_insert(cigar_op_new(h_right, 'H'), cigar_code->ops);
	} else {
	  h_right = 0;
	}

	//if (h_left > fq_read->length || h_left < 0) { exit(-1); }		

	//printf("FINAL H_LEFT = %i, H_RIGHT = %i\n", h_left, h_right);
	int len_read = fq_read->length - h_left - h_right;
	memcpy(query, &query_map[h_left], len_read);
	query[len_read] = '\0';

	memcpy(quality, &fq_read->quality[h_left], len_read);
	quality[len_read] = '\0';

	//printf("* * * %s * * *\n", fq_read->id);

	//int header_len = strlen(fq_read->id); 
	//char header_id[header_len + 1];
	//get_to_first_blank(fq_read->id, header_len, header_id);
	//char *header_match = (char *)malloc(sizeof(char)*header_len);	
	//memcpy(header_match, header_id, header_len);

	//printf("SCORE: %i\n", );

	if (!cigar_code_validate_(fq_read, cigar_code)) {
	  //meta_alignment_complete_free(meta_alignment);
	  //free(header_match);
	  continue;
	}

	alignment = alignment_new();
	alignment_init_single_end(strdup(fq_read->id),
				  strdup(query),
				  strdup(quality),
				  first_cal->strand, first_cal->chromosome_id - 1, start_mapping - 1,
				  strdup(new_cigar_code_string(cigar_code)),//strdup(cigar_fake)
				  cigar_code_get_num_ops(cigar_code),//1
				  cigar_code_score(cigar_code, fq_read->length),
				  1, 
				  (array_list_size(meta_alignments_list[i]) < 1),
				  cigar_code->distance, NULL, alignment);
	//alignment_print(alignment);
	
	//printf("Report : %s\n", new_cigar_code_string(cigar_code));

	array_list_insert(alignment, mapping_batch->mapping_lists[i]);
	meta_alignment_complete_free(meta_alignment);
	
      }
    } else {
      //data_type[i] = META_ALIGNMENT_TYPE;
      int num_items = array_list_size(meta_alignments_list[i]);
      int num_oks = 0;

      for (int m = num_items - 1; m >= 0; m--) {
	meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);	
	cigar_code_t *cigar_code = meta_alignment->cigar_code;
	int ok = 1;
	if (cigar_code_get_num_ops(cigar_code) <= 0) { 
	  ok = 0; 
	} else {
	  for (int j = 0; j < cigar_code->ops->size; j++) {
	    cigar_op_t *op = array_list_get(j, cigar_code->ops);
	    if (op->number <= 0)  { ok = 0; break; }
	  }
	}

	if (!ok) { 
	  array_list_remove_at(m, meta_alignments_list[i]); 
	  meta_alignment_complete_free(meta_alignment);
	}
      }

      if (array_list_size(meta_alignments_list[i]) > 0) {
	post_process_reads[i] = 1; 
	data_type[i] = META_ALIGNMENT_TYPE;
	for (int m = 0; m < array_list_size(meta_alignments_list[i]); m++) {
	  meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);
	  array_list_insert(meta_alignment, mapping_batch->mapping_lists[i]);
	}

	//pthread_mutex_lock(&mutex_sp);
	//extern size_t w3_r;
	//w3_r++;
	//pthread_mutex_unlock(&mutex_sp);

	//array_list_clear(meta_alignments_list[i], NULL);
	//pthread_mutex_lock(&mutex_sp);
	//insert_file_item_2(fq_read, meta_alignments_list[i], f_hc);
	//pthread_mutex_unlock(&mutex_sp);

	//pthread_mutex_lock(&mutex_sp);
	//extern size_t w2_r;
	//w2_r++;
	//file_write_items(fq_read, meta_alignments_list[i],
	//		 META_ALIGNMENT_TYPE, f_sa, f_hc, 0);
	//pthread_mutex_unlock(&mutex_sp);
	//array_list_clear(meta_alignments_list[i], (void*)meta_alignment_complete_free);
      } else {
	//printf(" ******************** No meta-alignment *****************\n");
	array_list_clear(mapping_batch->mapping_lists[i], (void*)NULL);
      }

    }

    if (data_type[i] == META_ALIGNMENT_TYPE) { continue; }

    size_t n_alignments = array_list_size(mapping_batch->mapping_lists[i]);
    int final_distance;
    for (size_t a = 0; a < n_alignments; a++) {
      alignment_t *alignment = (alignment_t *)array_list_get(a, mapping_batch->mapping_lists[i]);
      
      // set optional fields                                                             	
      optional_fields_length = 100;
      optional_fields = (char *) calloc(optional_fields_length, sizeof(char));

      final_distance = alignment->optional_fields_length;

      //float score_tmp = len_read * 0.5 - final_distance * 0.4;
      //int score_map = (int)(score_tmp * 254) / (len_read * 0.5);      
      //if (score_map < 0) { score_map = 0; }
      //else if (score_map > 254) { score_map = 254; }      
      p = optional_fields;
      AS = alignment->map_quality;
      alignment->map_quality = 255;

      sprintf(p, "ASi");
      p += 3;
      memcpy(p, &AS, sizeof(int));
      p += sizeof(int);
      
      sprintf(p, "NHi");
      p += 3;
      memcpy(p, &n_alignments, sizeof(int));
      p += sizeof(int);      

      sprintf(p, "NMi");
      p += 3;
      memcpy(p, &final_distance, sizeof(int));
      p += sizeof(int);	
      
      alignment->optional_fields_length = optional_fields_length;
      alignment->optional_fields = optional_fields;
    }

  }

  //printf("WK_1ph: =============== REPORT ALIGNMENTS END =====================\n"); 
  array_list_t *new_fq_batch = array_list_new(num_reads, 
					      1.25f, 
					      COLLECTION_MODE_ASYNCHRONIZED);
  int num_new_reads = 0;
 
  if (pair_mode == PAIRED_END_MODE || 
      pair_mode == MATE_PAIR_MODE) {
    //Pair end or Mate pair mode
    for (i = 0; i < num_reads; i +=2) {
      fastq_read_t *fq_read0 = array_list_get(i, mapping_batch->fq_batch);
      fastq_read_t *fq_read1 = array_list_get(i + 1, mapping_batch->fq_batch);
      if (data_type[i]     == CAL_TYPE || 
	  data_type[i + 1] == CAL_TYPE || 
	  data_type[i]     == META_ALIGNMENT_TYPE || 
	  data_type[i + 1] == META_ALIGNMENT_TYPE) {
	
	int mode;
	if (data_type[i] == CAL_TYPE || 
	    data_type[i + 1] == CAL_TYPE) {
	  mode = 0;
	} else {
	  mode = 1;
	}

	pthread_mutex_lock(&mutex_sp);
	file_write_items(fq_read0, mapping_batch->mapping_lists[i],
			 data_type[i], f_sa, f_hc, mode);
	file_write_items(fq_read1, mapping_batch->mapping_lists[i + 1],
			 data_type[i + 1], f_sa, f_hc, mode);
	pthread_mutex_unlock(&mutex_sp);

	if (data_type[i] == CAL_TYPE) {
	  array_list_clear(mapping_batch->mapping_lists[i], (void *)cal_free);
	} else if (data_type[i] == META_ALIGNMENT_TYPE) {
	  array_list_clear(mapping_batch->mapping_lists[i], (void *)meta_alignment_complete_free);
	} else {
	  array_list_clear(mapping_batch->mapping_lists[i], (void *)alignment_free);
	}
     
	if (data_type[i + 1] == CAL_TYPE) {
	  array_list_clear(mapping_batch->mapping_lists[i + 1], (void *)cal_free);
	} else if (data_type[i + 1] == META_ALIGNMENT_TYPE) {
	  array_list_clear(mapping_batch->mapping_lists[i + 1], (void *)meta_alignment_complete_free);
	} else {
	  array_list_clear(mapping_batch->mapping_lists[i + 1], (void *)alignment_free);
	}
      }

      array_list_free(meta_alignments_list[i], NULL);
      //if (rev_comp[i]) { free(rev_comp[i]); }
      free(scores_ranking[i]);

      array_list_free(meta_alignments_list[i + 1], NULL);
      //if (rev_comp[i + 1]) { free(rev_comp[i + 1]); }
      free(scores_ranking[i + 1]);

      if (post_process_reads[i] == 1 || post_process_reads[i + 1] == 1) {
	//pthread_mutex_lock(&mutex_sp);
	//extern size_t tot_reads_out;
	//tot_reads_out += 2;
	//pthread_mutex_unlock(&mutex_sp);

	fastq_read_free(fq_read0);
	array_list_free(mapping_batch->mapping_lists[i], NULL);

	fastq_read_free(fq_read1);
	array_list_free(mapping_batch->mapping_lists[i + 1], NULL);
      } else {
	//pthread_mutex_lock(&mutex_sp);
	//extern size_t tot_reads_in;
	//tot_reads_in += 2;
	//pthread_mutex_unlock(&mutex_sp);

	array_list_set_flag(1, mapping_batch->mapping_lists[i]);
	mapping_batch->mapping_lists[num_new_reads++] = mapping_batch->mapping_lists[i];
	array_list_insert(fq_read0, new_fq_batch);

	array_list_set_flag(1, mapping_batch->mapping_lists[i + 1]);
	mapping_batch->mapping_lists[num_new_reads++] = mapping_batch->mapping_lists[i + 1];
	array_list_insert(fq_read1, new_fq_batch);
      }
    }
  } else {
    //Single end mode
    for (i = 0; i < num_reads; i++) {
      array_list_free(meta_alignments_list[i], NULL);
      //if (rev_comp[i]) { free(rev_comp[i]); }
      free(scores_ranking[i]);

      fq_read = array_list_get(i, mapping_batch->fq_batch);      
      if (data_type[i] == CAL_TYPE || 
	  data_type[i] == META_ALIGNMENT_TYPE) {
	int mode;
	if (data_type[i] == CAL_TYPE) {
	  mode = 0;
	} else {
	  mode = 1;
	}
	pthread_mutex_lock(&mutex_sp);
	file_write_items(fq_read, mapping_batch->mapping_lists[i],
			 data_type[i], f_sa, f_hc, mode);
	pthread_mutex_unlock(&mutex_sp);
	
	if (data_type[i] == CAL_TYPE) {
	  array_list_clear(mapping_batch->mapping_lists[i], (void *)cal_free);
	} else if (data_type[i] == META_ALIGNMENT_TYPE) {
	  array_list_clear(mapping_batch->mapping_lists[i], (void *)meta_alignment_complete_free);
	}
      }
      //array_list_clear(alignments_list, alignment_free);      
      if (post_process_reads[i] == 0) {
	mapping_batch->mapping_lists[num_new_reads++] =  mapping_batch->mapping_lists[i];
	array_list_insert(fq_read, new_fq_batch);
      } else {
	fastq_read_free(fq_read);
	array_list_free(mapping_batch->mapping_lists[i], NULL);
      }      
    }
  }

  free(scores_ranking);
  free(post_process_reads);
  free(data_type);
  array_list_free(mapping_batch->fq_batch, NULL);
  mapping_batch->fq_batch = new_fq_batch;

  array_list_free(seeds_list, NULL);
  array_list_free(final_positions, NULL);
  free(new_targets);
  sw_multi_output_free(output);
  array_list_free(cals_targets, NULL);

  
  free(meta_alignments_list);
  //free(rev_comp);

  LOG_DEBUG("========= SPLICE JUNCTION SEARCH END =========\n");

  return RNA_POST_PAIR_STAGE;

}

int apply_rna_last(sw_server_input_t* input_p, batch_t *batch) {

  int min_intron_size = input_p->min_intron_size;
  size_t max_intron_size = input_p->max_intron_size;
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  sw_optarg_t *sw_optarg = &input_p->sw_optarg;
  genome_t *genome = input_p->genome_p;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  size_t num_targets = mapping_batch->num_targets;
  metaexons_t *metaexons = input_p->metaexons;
  cal_optarg_t *cal_optarg = input_p->cal_optarg_p;
  avls_list_t *avls_list = input_p->avls_list;
  bwt_optarg_t *bwt_optarg = input_p->bwt_optarg_p;
  bwt_index_t *bwt_index = input_p->bwt_index_p;
  linked_list_t *buffer = input_p->buffer;
  linked_list_t *buffer_hc = input_p->buffer_hc;
  int pair_mode = input_p->pair_mode;
  FILE *f_sa = input_p->f_sa;
  FILE *f_hc = input_p->f_hc;
  //fprintf(stderr, "APPLY RNA LAST START... %i\n", num_reads);

  array_list_t *cals_list, *fusion_cals, *fusion_cals_aux;
  cal_t *cal, *cal_prev, *cal_next, *first_cal, *last_cal;
  fastq_read_t *fq_read;
  linked_list_iterator_t itr;  
  seed_region_t *s, *s_prev, *s_next;
  cigar_code_t *cigar_code, *cigar_code_prev, *cigar_code_aux;
  cigar_code_t *alig_cigar_code;
  cigar_op_t *cigar_op_start, *cigar_op_end, *cigar_op, *cigar_op_prev, *cigar_op_aux;
  int *new_targets = (int *)calloc(num_reads, sizeof(int));
  array_list_t *merge_cals;
  linked_list_t *linked_list;
  seed_region_t *seed_region;
  //float *cals_score = (float *)calloc(100, sizeof(float));
  float score;
  char reference[2048];
  //char reference_prev[2048];
  //char reference_next[2048];
  //char reference_aux[2048];
  char query[2048];
  //char query_revcomp[2048];
  alignment_t *alignment;
  char q[2048];
  char r[2048];

  //char **rev_comp = (char **)calloc(num_reads, sizeof(char *));

  //fusion_coords_t *extrem_coords[2*40*num_reads];
  //fusion_coords_t *sp_coords[2*40*num_reads];
  //cigar_code_t *extrem_cigars[2*40*num_reads];
  //cigar_code_t *sp_cigars[2*40*num_reads];

  char *sequence;
  char *query_ref;
  char *quality_map, *query_map;
  //float scores_ranking[mapping_batch->num_allocated_targets][50];
  float *cals_score;
  char cigar_str[1024];

  cigar_op_t *first_op;
  char *match_seq, *match_qual, *optional_fields, *p;
  int match_start, match_len, optional_fields_length, AS;
  float norm_score;

  int delete_not_cigar;
  int padding_left, padding_right, len_query;  
  int num_sw = 0;
  int num_sw_sp = 0;
  int num_sw_ext = 0;
  size_t num_new_targets = 0;
  int coverage;
  int read_length;
  int cal_pos = 0;
  int number_of_best;
  //int num_extrem_ok = 0;
  //int num_sp_ok = 0;

  //Change to report most alignments
  int seed_err_size = 20;
  int n_alignments = 1;
  int target;
  int c, exact_nt;
  size_t genome_start, genome_end, genome_start2, genome_end2, read_start, read_end;
  int flank = 20;
  int n_delete;
  int sw_pos;
  int seeds_nt; 
  int e1;
  int e2;
  int ref_pos;
  int nt_2;
  
  int start_seeding, end_seeding;
  int lim_start, lim_end;

  //int min_intron_size = 40;

  metaexon_t *first_metaexon, *last_metaexon;
  int l_found, r_found; 
  size_t num_cals;
  register size_t t;
  register int i, j;

  int meta_type;
  meta_alignment_t *meta_alignment;
  sw_item_t *sw_item;  
  sw_depth_t sw_depth;
  sw_depth.depth = 0;
  sw_multi_output_t *output = sw_multi_output_new(MAX_DEPTH);
  //array_list_t **sw_items_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  //array_list_t **alignments_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  array_list_t **meta_alignments_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  array_list_t *seeds_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *final_positions = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED); 
  int make_seeds = 0;

  array_list_t *cals_targets = array_list_new(50, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  int *post_process_reads = (int *)calloc(num_reads, sizeof(int));
  //float scores_ranking[num_reads][100];
  float **scores_ranking = (float **)calloc(num_reads, sizeof(float *));//[num_reads][200];
  int read_nt;
  int *data_type  = (int *)calloc(num_reads, sizeof(int));
  struct timeval t_start, t_end;
  double time_s = 0;

  int from_single_anchors;
  int seed_size = cal_optarg->seed_size;

  //Delete!!!!
  //array_list_free(cals_targets, NULL);
  //printf("workflow w2 %i\n", w2_r);
  //  exit(-1);
  //1st - MAP THE EASY READS  

  for (i = 0; i < num_reads; i++) {
    meta_alignments_list[i] = array_list_new(20, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    fq_read = array_list_get(i, mapping_batch->fq_batch);
    cals_list = mapping_batch->mapping_lists[i];
    scores_ranking[i] = (float *)calloc(200, sizeof(float));
    from_single_anchors = 0;

    //printf("WK_2ph-Process: == (%i CALs)Read %s ==\n", array_list_size(cals_list), fq_read->id);
    
    //array_list_clear(mapping_batch->mapping_lists[i], (void*)cal_free);
    //continue;

    /* 
    if (array_list_size(cals_list) == 0) { continue; }
    if (array_list_get_flag(cals_list) == BITEM_SINGLE_ANCHORS) {
      array_list_clear(cals_list, cal_free);
    } else if (array_list_get_flag(cals_list) == BITEM_CALS) {
      for (int t = 0; t < array_list_size(cals_list); t++) {
	array_list_t *fusion_list = array_list_get(t, cals_list);
	array_list_free(fusion_list, cal_free);
      }
      array_list_clear(cals_list, NULL);
    } else if (array_list_get_flag(cals_list) == BITEM_META_ALIGNMENTS) {
      array_list_clear(cals_list, meta_alignment_complete_free);
    }
    continue;
    */

    if (array_list_get_flag(cals_list) != BITEM_SINGLE_ANCHORS) { 
      if (array_list_get_flag(cals_list) == BITEM_META_ALIGNMENTS) {
	post_process_reads[i] = 1; 
	data_type[i] = META_ALIGNMENT_TYPE;
      }
      continue;
    }
    
    //if (array_list_size(cals_list) > 0) {
    assert(array_list_size(cals_list) > 0);

    //if (array_list_get_flag(cals_list) == BITEM_SINGLE_ANCHORS) {
      //SINGLE ANCHORS FOUND

    int map_read = 0, map = 0;      
    int read_nt;
    int seeds_process = 0;
    cigar_code_t *cigar_code_res;

    for (int p = 0; p < array_list_size(cals_list); p++) {
      cal = array_list_get(p, cals_list);
      //cal_print(cal);
      cigar_code_res = NULL;
      make_seeds = 0;
      if (cal->strand == 1) {
	//if (!rev_comp[i]) {
	//rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	//strcpy(rev_comp[i], fq_read->sequence);
	//seq_reverse_complementary(rev_comp[i], fq_read->length);
	//}
	//query_map = rev_comp[i];
	query_map = fq_read->revcomp;
      } else {
	query_map = fq_read->sequence;
      }

      s_prev = linked_list_get_last(cal->sr_list);
      if (s_prev->read_start == 0) {
	read_nt = fq_read->length - (s_prev->read_end + 1);
      } else {
	read_nt = s_prev->read_start;
      }
      if (s_prev->read_start == 0) {
	cigar_code_res = fill_extrem_gap(query_map, 
					 cal,
					 FILL_GAP_RIGHT,
					 genome,
					 metaexons, 
					 avls_list);
      } else {
	cigar_code_res = fill_extrem_gap(query_map, 
					 cal,
					 FILL_GAP_LEFT,
					 genome,
					 metaexons, 
					 avls_list);	
      }
      if (cigar_code_res) {
	meta_alignment = meta_alignment_new();
	meta_alignment_insert_cal(cal, meta_alignment);
	  
	if (s_prev->read_start == 0) {
	  meta_alignment->cigar_right = cigar_code_res;
	} else {
	  meta_alignment->cigar_left = cigar_code_res;
	}
	    
	meta_alignment_fill_gaps(META_ALIGNMENT_RIGHT,
				 meta_alignment, query_map, genome,
				 sw_optarg, output, metaexons, 
				 &sw_depth, avls_list, min_intron_size);
	    
	array_list_insert((void *)meta_alignment, meta_alignments_list[i]);
	    
      }
    }

    //if (array_list_size(meta_alignments_list[i]) <= 0) {	
    //array_list_clear(mapping_batch->mapping_lists[i], (void*)cal_free);
    //} 
 
    //continue;

    if (array_list_size(meta_alignments_list[i]) <= 0) {	
      //Make Seeding and Caling
      array_list_clear(mapping_batch->mapping_lists[i], (void*)cal_free);
      //array_list_clear(cals_list, NULL);
      array_list_t *new_cals_list = mapping_batch->mapping_lists[i];
      //int seed_size = 16;
      num_cals = bwt_generate_cals(fq_read->sequence, 
				   seed_size, bwt_optarg,
				   cal_optarg,
				   bwt_index, new_cals_list, 
				   genome->num_chromosomes + 1);
      //filter-incoherent CALs
      int founds[num_cals], found = 0;
      for (size_t j = 0; j < num_cals; j++) {
	founds[j] = 0;
	cal = array_list_get(j, new_cals_list);
	LOG_DEBUG_F("\tcal %i of %i: sr_list size = %i (cal->num_seeds = %i) %i:%lu-%lu\n", 
		    j, num_cals, cal->sr_list->size, cal->num_seeds,
		    cal->chromosome_id, cal->start, cal->end);
	if (cal->sr_list->size > 0) {
	  int start = 0;
	  size_t genome_start = 0;
	  int  first = 1;
	  for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
	    seed_region_t *s = list_item->item;
	      
	    LOG_DEBUG_F("\t\t:: star %lu > %lu s->read_start\n", start, s->read_start);
	    if (start > s->read_start || 
		s->read_start >= s->read_end) {
	      LOG_DEBUG("\t\t\t:: remove\n");
	      found++;
	      founds[j] = 1;
	    }
	    if (!first && 
		((s->genome_start < genome_start) || 
		 (s->genome_start - genome_start) > 2*fq_read->length)) {
	      //printf("Remove (genome_start = %i s->genome_start = %i)\n", genome_start, s->genome_start);
	      //cal_print(cal);
	      found++;
	      founds[j] = 1;
	    }

	    first = 0;
	    start = s->read_end + 1;
	    genome_start = s->genome_end + 1;
	  }
	} else {
	  found++;
	  founds[j] = 1;
	}
      }

      array_list_t *cal_list_aux;
      if (found) {
	int min_seeds = 100000;
	int max_seeds = 0;
	cal_list_aux = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
	for (size_t j = 0; j < num_cals; j++) {
	  if (!founds[j]) {
	    cal = array_list_get(j, new_cals_list);
	    cal->num_seeds = cal->sr_list->size;
	    if (cal->num_seeds > max_seeds) max_seeds = cal->num_seeds;
	    if (cal->num_seeds < min_seeds) min_seeds = cal->num_seeds;
	    array_list_insert(cal, cal_list_aux);
	    array_list_set(j, NULL, new_cals_list);
	  }
	}

	array_list_free(new_cals_list, (void *) cal_free);
	num_cals = array_list_size(cal_list_aux);
	new_cals_list = cal_list_aux;
	mapping_batch->mapping_lists[i] = cal_list_aux;
      }

      if (num_cals) {
	int max = 100;
	if (num_cals > max) {
	  int select_cals = num_cals - max;
	  for(int j = num_cals - 1; j >= max; j--) {
	    cal_free(array_list_remove_at(j, new_cals_list));
	  }
	}
      }
    

      cals_list = mapping_batch->mapping_lists[i];
      
      if (array_list_size(cals_list) > 0) { 
	//if (array_list_get_flag(cals_list) == BITEM_CALS) {    
	//CALS FOUND
	//printf("\tWK_2ph: -- CALS PROCESS -- %i\n", array_list_size(cals_list));     
	//if (array_list_size(meta_alignments_list[i]) > 0) {	
	//array_list_clear(mapping_batch->mapping_lists[i], (void*)NULL);
	//}// else {
	//if (!from_single_anchors) {
	//cals_targets = cals_list;
	array_list_clear(cals_targets, (void*)NULL);
	merge_and_filter_cals(cals_targets, cals_list, fq_read, 
			      bwt_optarg, bwt_index, genome, 
			      scores_ranking[i]);	    

	//}	
	//start_timer(t_start); 
	//===== Step-2: Generate CALs Score =====//
	LOG_DEBUG("STEP-2: GENERATE CALS SCORE");
	for (int s = 0; s < array_list_size(cals_targets); s++) {
	  fusion_cals = array_list_get(s, cals_targets);
	  scores_ranking[i][s] = generate_cals_score(fusion_cals, fq_read->length);
	  
	  metaexon_t *first_metaexon;
	  cal = array_list_get(0, fusion_cals);
	}
	//===== Step-2: END =====//
	
	size_t number_of_best = array_list_size(cals_targets);
	
	//===== Step-3: Ranking CALs by score (the n best)=====//
	LOG_DEBUG("STEP-3: ORDER CALS");
	float score;
	for (int s = 0; s < number_of_best - 1; s++) {
	  for (int s1 = s + 1; s1 < array_list_size(cals_targets); s1++) {
	    if (scores_ranking[i][s1] > scores_ranking[i][s]) {
	      array_list_swap(s, s1, cals_targets);
	      score = scores_ranking[i][s1];
	      scores_ranking[i][s1] = scores_ranking[i][s];
	      scores_ranking[i][s] = score;
	    }
	  }
	}
	//===== Step-3: END =====//

	//int MAX_CALS_PROCESS = 10;
	//size_t seed_size = 16;
	float best_score = scores_ranking[i][0];
	int limit;	

	if (array_list_size(cals_targets) > 5) {
	  for (int j = array_list_size(cals_targets) - 1; j >= 5; j--) { 
	    array_list_t *fusion_list = array_list_remove_at(j, cals_targets);
	    array_list_free(fusion_list, (void *)cal_free);
	  }
	}
	
	
	limit = array_list_size(cals_targets);// > 5 ? 5 : array_list_size(cals_targets);
	
	int num_process = 0;
	//fprintf(stderr, "%i vs %i : %s\n", limit, array_list_size(cals_targets), 
	//	fq_read->id);
	for (int j = 0; j < limit; j++) { 
	  meta_type = -1;
	  array_list_t *fusion_list = array_list_get(j, cals_targets);
		  
	  if (scores_ranking[i][j] == 0) { 	    
	    array_list_free(fusion_list, (void *)cal_free);
	    continue; 
	  }
	  /*
	  printf("==== FUSION CALS PROCESS ====\n");
	  for (int t = 0; t < array_list_size(fusion_list); t++) {
	    cal_t *cal_aux = array_list_get(t, fusion_list);
	    cal_print(cal_aux);
	  }
	  printf("==== ------------------- ====\n");
	  */
	  first_cal = array_list_get(0, fusion_list);
	  last_cal = array_list_get(array_list_size(fusion_list) - 1, fusion_list);
	  if (first_cal->strand == 1) {
	    //if (!rev_comp[i]) {
	    //rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	    //strcpy(rev_comp[i], fq_read->sequence);
	    //seq_reverse_complementary(rev_comp[i], fq_read->length);
	    //}
	    //query_map = rev_comp[i];
	    query_map = fq_read->revcomp;
	  } else {
	    query_map = fq_read->sequence;
	  }

	  s_prev = linked_list_get_first(first_cal->sr_list);    
	  s_next = linked_list_get_last(last_cal->sr_list);		    	  

	  num_process++;
	  //printf("==== ------------------- ====\n");
	  if (array_list_size(fusion_list) > 1) { 
	    //printf("FUSION CALS REPORT\n");
	    cal_t *first_cal = array_list_get(0, fusion_list);
	    cal_t *last_cal  = array_list_get(array_list_size(fusion_list) - 1, fusion_list);
	   
	    //TODO: IF WE HAVE MORE THAN TWO CALS SEARCH SINGLE ANCHOR
	    cigar_code = search_double_anchors_cal(query_map,
						   first_cal, last_cal,
						   metaexons, genome,
						   fq_read, &meta_type, avls_list);
	    if (cigar_code != NULL) {
	      //printf("FOUND! %s\n", new_cigar_code_string(cigar_code));
	      meta_alignment = meta_alignment_new();
	      meta_alignment_insert_cal(first_cal, meta_alignment);
	      if (meta_type == META_ALIGNMENT_LEFT) {
		meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
	      } else {
		meta_alignment_insert_cal(last_cal, meta_alignment);
		meta_alignment_insert_cigar(cigar_code, CIGAR_SIMPLE_MIDDLE, 0, meta_alignment);
	      }

	      meta_alignment_fill_gaps(meta_type,
				       meta_alignment, query_map, genome,
				       sw_optarg, output, metaexons, 
				       &sw_depth, avls_list, min_intron_size);	
	    } else {
	      //printf("NOT FOUND\n");
	      meta_alignment = meta_alignment_cals_new(fusion_list);  
	      for (int c = 1; c < array_list_size(fusion_list); c++) { 
		last_cal = array_list_get(c, fusion_list);
		avl_node_t *node_avl_start, *node_avl_end;
		seed_region_t *s_prev = linked_list_get_last(first_cal->sr_list);    
		seed_region_t *s_next = linked_list_get_first(last_cal->sr_list);	
		int distance_aux;
		size_t sp_start, sp_end;
		int sp_type;	
		int nt = search_simple_splice_junction(s_prev, s_next, 
						       last_cal->chromosome_id, 
						       last_cal->strand,
						       query_map, genome, 
						       &sp_start, &sp_end,
						       &sp_type,
						       &distance_aux);

		if (nt) {
		  cigar_code = cigar_code_new();
		  cigar_code->distance = distance_aux;
		  cigar_code_append_new_op(nt, 'N', cigar_code);
		
		  seed_region_t *s_prev_aux = linked_list_get_last(first_cal->sr_list);    
		  seed_region_t *s_next_aux = linked_list_get_first(last_cal->sr_list);

		  int sp_strand = (sp_type == CT_AC_SPLICE ? 1 : 0);

		  allocate_start_node(first_cal->chromosome_id - 1,
				      sp_strand,
				      sp_start,
				      sp_end,
				      sp_start,
				      sp_end,
				      FROM_READ,
				      sp_type,
				      NULL, 
				      &node_avl_start,
				      &node_avl_end, 
				      avls_list);		  

		  assert(s_prev_aux->genome_start < node_avl_start->position);
		  assert(node_avl_end->position < s_next_aux->genome_end);
		  
		  metaexon_insert(first_cal->strand, first_cal->chromosome_id - 1,
				  s_prev_aux->genome_start, node_avl_start->position, 40,
				  METAEXON_RIGHT_END, node_avl_start,
				  metaexons);
		
		  metaexon_insert(last_cal->strand, last_cal->chromosome_id - 1,
				  node_avl_end->position, s_next_aux->genome_end, 40,
				  METAEXON_LEFT_END, node_avl_end,
				  metaexons);
		  
		  meta_alignment_insert_cigar(cigar_code, CIGAR_SIMPLE_MIDDLE, c - 1, meta_alignment);
		
		} else {
		  
		  info_sp_t *info_sp = sw_reference_splice_junction(first_cal, last_cal,
								    query_map, genome,
								    q, r);
		  //New sw item. Storing data...
		  sw_item = sw_item_new(SP_SW, i, 0, c - 1,
					first_cal, last_cal, 
					meta_alignment, NULL, 
					NULL, info_sp);	      
		  //Insert item... and process if depth is full
		  sw_depth_insert(q, r, sw_item,
				  sw_optarg, output,
				  avls_list, metaexons,
				  &sw_depth, genome, min_intron_size); 

		  //meta_alignment_insert_cigar(NULL, CIGAR_SIMPLE_MIDDLE, c - 1, meta_alignment);
		}

		first_cal = last_cal;

	      }

	      //array_list_free(fusion_cals, NULL);
	      meta_alignment_fill_gaps(META_ALIGNMENT_MIDDLE,
				       meta_alignment, query_map, genome,
				       sw_optarg, output, metaexons, 
				       &sw_depth, avls_list, min_intron_size);	 
	    }

	  } else {
	    //printf(":::: SINGLE MAP (%i)::::\n", array_list_size(cals_list));
	    cal = array_list_get(0, fusion_list);
	    //cal_print(cal);
	  
	    meta_alignment = meta_alignment_new();
	    meta_alignment_insert_cal(cal, meta_alignment);
	    meta_alignment_fill_gaps(META_ALIGNMENT_MIDDLE,
				     meta_alignment, query_map, genome,
				     sw_optarg, output, metaexons, 
				     &sw_depth, avls_list, min_intron_size); 
	  }
	
	  array_list_insert(meta_alignment, meta_alignments_list[i]); 
	

	  cal_t *first_cal = array_list_get(0, fusion_list);
	  cal_t *last_cal  = array_list_get(array_list_size(fusion_list) - 1, fusion_list);	    
	  seed_region_t *s_prev = linked_list_get_first(first_cal->sr_list);
	  seed_region_t *s_next = linked_list_get_last(last_cal->sr_list);
	  cigar_code_t *cc_left, *cc_right;
	  
	  if (s_prev->read_start != 0) {
	    cc_left = fill_extrem_gap(query_map, 
				      first_cal,
				      FILL_GAP_LEFT,
				      genome,
				      metaexons,
				      avls_list); 
	    meta_alignment_insert_cigar(cc_left, CIGAR_ANCHOR_RIGHT, 0, meta_alignment);
	  }
	
	  if (meta_type != META_ALIGNMENT_LEFT &&
	      s_next->read_end != fq_read->length - 1) {
	    cc_right = fill_extrem_gap(query_map, 
				       last_cal,
				       FILL_GAP_RIGHT,
				       genome,
				       metaexons, 
				       avls_list); 
	    meta_alignment_insert_cigar(cc_right, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
	  }
	  array_list_free(fusion_list, NULL);	  
	}
      }
    }
  
    array_list_clear(mapping_batch->mapping_lists[i], (void*)NULL);

    //printf("::::AFTER CALS %i\n", num_process);
    //stop_timer(t_start, t_end, time_s);
    //extern double seeding_time_2;
    //extern pthread_mutex_t mutex_sp; 
    //pthread_mutex_lock(&mutex_sp);
    //seeding_time_2 += time_s;
    //time_s = 0;
    //pthread_mutex_unlock(&mutex_sp);
    
    
    /*else if (array_list_get_flag(cals_list) == BITEM_META_ALIGNMENTS) {
      
    //META ALIGNMENTS
    //printf("\tWK_2ph: -- META-ALIGNMENTS PROCESS -- \n");
    //fprintf(stderr, "\tMETA ALIGNMENTS\n");
    for (int j = 0; j < array_list_size(cals_list); j++) {
    meta_alignment = array_list_get(j, cals_list);
    
    //printf("SCORE IN : %i\n", meta_alignment->score);
    array_list_t *fusion_list = meta_alignment->cals_list;
    cal_t *first_cal = array_list_get(0, fusion_list);
    cal_t *last_cal  = array_list_get(array_list_size(fusion_list) - 1, fusion_list);	    
    seed_region_t *s_prev = linked_list_get_first(first_cal->sr_list);
    seed_region_t *s_next = linked_list_get_last(last_cal->sr_list);
    cigar_code_t *cc_left, *cc_right;
    
    if (first_cal->strand == 1) {
    if (!rev_comp[i]) {
    rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
    strcpy(rev_comp[i], fq_read->sequence);
    seq_reverse_complementary(rev_comp[i], fq_read->length);
    }
    query_map = rev_comp[i];
    } else {
    query_map = fq_read->sequence;
    }
	  
    if (s_prev->read_start != 0 && 
    meta_alignment->cigar_left == NULL) {
    metaexon_search(first_cal->strand, first_cal->chromosome_id - 1,
    first_cal->start, first_cal->end, &first_metaexon,
    metaexons);
    cc_left = fill_extrem_gap(query_map, 
    first_cal,
    FILL_GAP_LEFT,
    genome,
    metaexons, 
    avls_list); 
    assert(meta_alignment != NULL);
    meta_alignment_insert_cigar(cc_left, CIGAR_ANCHOR_RIGHT, 0, meta_alignment);
    }
    
    if (s_next->read_end != fq_read->length - 1 &&
    meta_alignment->cigar_right == NULL) {
    metaexon_search(last_cal->strand, last_cal->chromosome_id - 1,
    last_cal->start, last_cal->end, &first_metaexon,
    metaexons);	    
    cc_right = fill_extrem_gap(query_map, 
    last_cal,
    FILL_GAP_RIGHT,
    genome,
    metaexons, 
    avls_list); 
    //printf("RESULT CIGAR %s\n", new_cigar_code_string(cc_right));
    meta_alignment_insert_cigar(cc_right, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
    }	
    }
    } else {
    //NO CALS FOUND
    //printf("\tWK_2ph: -- NOT CALS FOUND --\n");
    //fprintf(stderr, "\tNO CALS\n");
    }*/
    //}
    //array_list_set_flag(BITEM_SINGLE_ANCHORS, mapping_batch->mapping_lists[i]);
  }
  
  sw_depth_process(sw_optarg, output, 
		   &sw_depth, avls_list, metaexons, 
		   SW_FINAL, genome, min_intron_size);

  //fprintf(stderr, "CLOSE META ALIGNMENTS\n");
  for (int i = 0; i < num_reads; i++) {
    fq_read = array_list_get(i, mapping_batch->fq_batch);

    if (array_list_size(meta_alignments_list[i]) <= 0 ||
	array_list_get_flag(mapping_batch->mapping_lists[i]) != BITEM_SINGLE_ANCHORS) { 
      continue;
    }

    //printf(".SECOND. : %s\n", fq_read->id);
    //if (array_list_get_flag(mapping_batch->mapping_lists[i]) != BITEM_META_ALIGNMENTS) {
    //array_list_clear(mapping_batch->mapping_lists[i], (void*)NULL);    
    //printf("<<<<CLOSE META (%i) %s\n", array_list_size(meta_alignments_list[i]), fq_read->id);
    for (int m = 0; m < array_list_size(meta_alignments_list[i]); m++) {
      meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);
      //printf("Status %i == %i\n", meta_alignment_get_status(meta_alignment), META_OPEN);
      if (meta_alignment_get_status(meta_alignment) == META_OPEN ) {
	meta_alignment_close(meta_alignment);
	//printf("CIGAR CLOSE %i: %s\n", m, new_cigar_code_string(meta_alignment->cigar_code));
      }
    }

    //} else {      
    //for (int j  = array_list_size(mapping_batch->mapping_lists[i]) - 1; j >= 0; j--) {
    //meta_alignment_t *meta_alignment = array_list_remove_at(j, mapping_batch->mapping_lists[i]);
    //meta_alignment_close(meta_alignment);
    //printf("CLOSE META %i : %s\n", j, new_cigar_code_string(meta_alignment->cigar_code));
    //meta_alignment_calculate_score(meta_alignment);
    //printf("\n2: SCORE-META %i\n", meta_alignment->score);
    //array_list_insert(meta_alignment, meta_alignments_list[i]);
    //}

    int no_map = 1;
    for (int m = 0; m < array_list_size(meta_alignments_list[i]); m++) {
      meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);
      if (meta_alignment->score == fq_read->length) {
	no_map = 0;
	break;
      }
    }
    
    meta_alignment_t *meta_alignment = array_list_get(0, meta_alignments_list[i]);
    if (no_map) { 
      meta_alignment->flag = 1;
      //printf("\t.SECOND-NO-MAP. :%s\n", fq_read->id);
    } else {
      meta_alignment->flag = 0;
      //printf("\t.SECOND-MAP. :%s\n", fq_read->id);
    }

  }

  size_t start_mapping;
  //printf("WK_2ph: =============== REPORT ALIGNMENTS =====================\n");
  //fprintf(stderr, "REPORT META ALIGNMENTS\n");
  for (int i = 0; i < num_reads; i++) {
    fq_read = array_list_get(i, mapping_batch->fq_batch);
    //printf("WK_2ph-Meta_Report:  %s >>>>\n", fq_read->id);
    //printf("Meta alignments flag value if (%i == 0|| %i != %i)\n",
    //   array_list_size(meta_alignments_list[i]),
    //	   array_list_get_flag(mapping_batch->mapping_lists[i]),
    //	   BITEM_SINGLE_ANCHORS);

    if (array_list_size(meta_alignments_list[i]) <= 0 ||
	array_list_get_flag(mapping_batch->mapping_lists[i]) != BITEM_SINGLE_ANCHORS) { 
      //printf("Continue...\n");
      continue;
    }

    //printf("QUE: %s\nQUA: %s\n", fq_read->sequence, fq_read->quality);
    assert(fq_read->id != NULL);

    char query[2048];
    char quality[2048];
    int map = 0;

    //meta_alignment_t *meta_alignment = array_list_get(0, meta_alignments_list[i]);
    //if (meta_alignment->score == fq_read->length) {
    meta_alignment_t *meta_alignment = array_list_get(0, meta_alignments_list[i]);      
    if (meta_alignment->flag == 0) {
      data_type[i] = ALIGNMENT_TYPE;
      for (int m = 0; m < array_list_size(meta_alignments_list[i]); m++) {
	//printf("Report meta %i\n", m);
	meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);      
	if (meta_alignment_get_status(meta_alignment) == META_CLOSE) { 
	  //if (fq_read->id[1] == 'r') { printf("::--> %s\n", fq_read->id); exit(-1); }
	  optional_fields_length = 0;
	  optional_fields = NULL;
	  
	  //printf("= SHOWING CALS  =\n");
	  //cal_print(first_cal);
	  //printf("= SHOWING ENDS  =\n");
	  cal_t *first_cal = meta_alignment_get_first_cal(meta_alignment);
	  s_prev = linked_list_get_first(first_cal->sr_list);
	  int h_left = s_prev->read_start;
	    
	  cal_t *last_cal = meta_alignment_get_last_cal(meta_alignment);
	  s_next = linked_list_get_last(last_cal->sr_list);
	  int h_right = (fq_read->length - 1) - s_next->read_end;
	    
	  if (first_cal->strand == 1) {
	    //query_map = rev_comp[i];
	    query_map = fq_read->revcomp;
	  } else {
	    query_map = fq_read->sequence;
	  }
	    
	  cigar_code = meta_alignment->cigar_code;	
	  assert(cigar_code != NULL);
	  start_mapping = first_cal->start;
	  int dsp = 0;
	  
	  //fprintf(stderr, "Read length %i - \n", fq_read->length);
	  if (h_left > 0 && 
	      meta_alignment->cigar_left == NULL &&
	      meta_alignment->type != META_ALIGNMENT_LEFT &&
	      meta_alignment->type != META_ALIGNMENT_RIGHT) {
	    //fprintf(stderr, "H_LEFT ->  %i -> seed %i\n", h_left, s_prev->read_start);
	    array_list_insert_at(0, cigar_op_new(h_left, 'H'), cigar_code->ops);
	  } else {
	    h_left = 0;
	    if (meta_alignment->cigar_left != NULL ) {
	      cigar_code_t *cigar_code = meta_alignment->cigar_left;
	      for (int c = 0; c < cigar_code->ops->size; c++) {
		cigar_op_t *op = array_list_get(c, cigar_code->ops);
		if (op->name == 'M' ||
		    op->name == 'N' ||
		    op->name == 'D') {
		  dsp += op->number;
		}
	      }
	    } 
	  }

	  start_mapping -= dsp;
	  //printf(" new_start = %lu\n", start_mapping);
	    
	  if (h_right > 0 &&
	      meta_alignment->cigar_right == NULL &&
	      meta_alignment->type != META_ALIGNMENT_LEFT &&
	      meta_alignment->type != META_ALIGNMENT_RIGHT) {
	    //fprintf(stderr, "H_RIGHT ->  %i -> seed %i\n", h_right,  s_next->read_end);
	    array_list_insert(cigar_op_new(h_right, 'H'), cigar_code->ops);
	  } else {
	    h_right = 0;
	  }
	    
	  //printf("FINAL H_LEFT = %i, H_RIGHT = %i\n", h_left, h_right);
	    
	  if (h_left > fq_read->length || h_left < 0) { exit(-1); }
	  if (h_left + h_right >= fq_read->length) { continue; }


	  if (!cigar_code_validate_(fq_read, cigar_code)) {
	    meta_alignment_complete_free(meta_alignment);
	    continue;
	  }
	    
	  int len_read = fq_read->length - h_left - h_right;
	  //printf("(REAL %i): %s (%i)\n", len_read, query, s);
	  memcpy(query, &query_map[h_left], len_read);
	  query[len_read] = '\0';
	    
	  memcpy(quality, &fq_read->quality[h_left], len_read);
	  quality[len_read] = '\0';
	    
	  //printf("* * * %s * * *\n", fq_read->id);
	  //fprintf(stderr, "* * * (%i) M E T A    A L I G N M E N T    R E P O R T    %s* * *\n", strlen(query), new_cigar_code_string(cigar_code));
	  alignment = alignment_new();
	    
	  //int header_len = strlen(fq_read->id); 
	  //char header_id[header_len + 1];
	  //get_to_first_blank(fq_read->id, header_len, header_id);
	  //char *header_match = (char *)malloc(sizeof(char)*header_len);
	  //if (header_match == NULL) { exit(-1); }
	  //memcpy(header_match, header_id, header_len);

	  //float score_tmp = len_read * 0.5 - cigar_code->distance * 0.4;
	  //int score_map = (int)(score_tmp * 254) / (len_read * 0.5);
	    
	  /*if (!cigar_code_validate(len_read, cigar_code)) {
	    char cigar_fake[512];
	    sprintf(cigar_fake, "%iM", fq_read->length);
	    //fprintf(stderr, "WK_2ph: * * * M E T A    A L I G N M E N T    R E P O R T    F A K E  [%i:%lu]  %s* * *\n", 
	    //	   first_cal->chromosome_id, first_cal->start, new_cigar_code_string(cigar_code));
	    //fprintf(stderr, "@@@@@@@@@ :%s\n", fq_read->id);
	    //fprintf(stderr, "@FAKE :%s\n", fq_read->id);
	    alignment_init_single_end(header_match, 
				      strdup(fq_read->sequence),
				      strdup(fq_read->quality),
				      first_cal->strand, first_cal->chromosome_id - 1, start_mapping,
				      strdup(cigar_fake),
				      1,
				      norm_score * 254, 1, (array_list_size(meta_alignments_list[i]) >= 1),
				      optional_fields_length, optional_fields, 0, alignment);	  
	    //fprintf(stderr, "OK INSERT\n");
	    } else {*/
	    //fprintf(stderr, "META ALIGNMENT REPORT %i: %s\n", m, new_cigar_code_string(cigar_code));
	    //printf("WK_2ph: * * * M E T A    A L I G N M E N T    R E P O R T  [%i:%lu]  %s* * *\n", 
	    //	   first_cal->chromosome_id, first_cal->start, new_cigar_code_string(cigar_code));
	  alignment_init_single_end(strdup(fq_read->id),
				    strdup(query)/*match_seq*/,
				    strdup(quality)/*match_qual*/,
				    first_cal->strand, first_cal->chromosome_id - 1, start_mapping - 1,
				    strdup(new_cigar_code_string(cigar_code))/*strdup(cigar_fake)*/,
				    cigar_code_get_num_ops(cigar_code)/*1*/,
				    cigar_code_score(cigar_code, fq_read->length), 
				    1, (array_list_size(meta_alignments_list[i]) < 1),
				    cigar_code->distance, NULL, alignment);
	    //fprintf(stderr, "OK INSERT\n");
	    //alignment_print(alignment);
	    //}
	  array_list_insert(alignment, mapping_batch->mapping_lists[i]);
	  //printf("Insert ok!\n");
	  map = 1;

	  meta_alignment_complete_free(meta_alignment);
	}
      }
    } else {
      int num_items = array_list_size(meta_alignments_list[i]);
      int num_oks = 0;

      for (int m = num_items - 1; m >= 0; m--) {
	meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);	
	cigar_code_t *cigar_code = meta_alignment->cigar_code;
	int ok = 1;
	if (cigar_code_get_num_ops(cigar_code) <= 0) { 
	  ok = 0; 
	} else {
	  for (int j = 0; j < cigar_code->ops->size; j++) {
	    cigar_op_t *op = array_list_get(j, cigar_code->ops);
	    if (op->number <= 0)  { ok = 0; break; }
	  }
	}

	if (!ok) { 
	  array_list_remove_at(m, meta_alignments_list[i]); 
	  meta_alignment_complete_free(meta_alignment);
	}
      }

      if (array_list_size(meta_alignments_list[i]) > 0) {
	//pthread_mutex_lock(&mutex_sp);
	//insert_file_item_2(fq_read, meta_alignments_list[i], f_hc);
	//pthread_mutex_unlock(&mutex_sp);
	//pthread_mutex_lock(&mutex_sp);
	//file_write_items(fq_read, meta_alignments_list[i],
	//		 META_ALIGNMENT_TYPE, f_sa, f_hc, 0);
	//pthread_mutex_unlock(&mutex_sp);

	post_process_reads[i] = 1; 
	data_type[i] = META_ALIGNMENT_TYPE;
	for (int m1 = 0; m1 < array_list_size(meta_alignments_list[i]); m1++) {
	  meta_alignment_t *meta_alignment = array_list_get(m1, meta_alignments_list[i]);
	  array_list_insert(meta_alignment, mapping_batch->mapping_lists[i]);
	}	
	//array_list_clear(meta_alignments_list[i], (void*)meta_alignment_complete_free);
	//array_list_clear(mapping_batch->mapping_lists[i], (void*)NULL);
      } else {
	array_list_clear(mapping_batch->mapping_lists[i], (void*)NULL);
      }
    }

    if (data_type[i] == META_ALIGNMENT_TYPE) { continue; }

    size_t n_alignments = array_list_size(mapping_batch->mapping_lists[i]);
    int final_distance;
    for (size_t a = 0; a < n_alignments; a++) {
      alignment_t *alignment = (alignment_t *)array_list_get(a, mapping_batch->mapping_lists[i]);
      
      // set optional fields                                                             	
      optional_fields_length = 100;
      optional_fields = (char *) calloc(optional_fields_length, sizeof(char));
      
      p = optional_fields;
      AS = alignment->map_quality;
      alignment->map_quality = 255;
      
      sprintf(p, "ASi");
      p += 3;
      memcpy(p, &AS, sizeof(int));
      p += sizeof(int);
      
      sprintf(p, "NHi");
      p += 3;
      memcpy(p, &n_alignments, sizeof(int));
      p += sizeof(int);
      
      sprintf(p, "NMi");
      p += 3;
      final_distance = alignment->optional_fields_length;
      memcpy(p, &final_distance, sizeof(int));
      p += sizeof(int);	
      
      alignment->optional_fields_length = optional_fields_length;
      alignment->optional_fields = optional_fields;
    }
  }

  array_list_t *new_fq_batch = array_list_new(num_reads, 
					      1.25f, 
					      COLLECTION_MODE_ASYNCHRONIZED);
  int num_new_reads = 0;
  if (pair_mode == PAIRED_END_MODE || 
      pair_mode == MATE_PAIR_MODE) {
    //Pair end or Mate pair mode
    for (i = 0; i < num_reads; i +=2) {
      fastq_read_t *fq_read0 = array_list_get(i, mapping_batch->fq_batch);
      fastq_read_t *fq_read1 = array_list_get(i + 1, mapping_batch->fq_batch);
      if (data_type[i]     == CAL_TYPE || 
	  data_type[i + 1] == CAL_TYPE || 
	  data_type[i]     == META_ALIGNMENT_TYPE || 
	  data_type[i + 1] == META_ALIGNMENT_TYPE) {
	
	int mode;
	if (data_type[i] == CAL_TYPE || 
	    data_type[i + 1] == CAL_TYPE) {
	  mode = 0;
	} else {
	  mode = 1;
	}

	pthread_mutex_lock(&mutex_sp);
	//printf("\n====\n(0)%s\n(1)%s\n====\n(0)%i-Insert num items %i\n(1)%i-Insert num items %i\n\n", fq_read0->id, fq_read1->id, data_type[i], array_list_size(mapping_batch->mapping_lists[i]), data_type[i+1], array_list_size(mapping_batch->mapping_lists[i + 1]));
	
	file_write_items(fq_read0, mapping_batch->mapping_lists[i],
			 data_type[i], f_sa, f_hc, mode);
	file_write_items(fq_read1, mapping_batch->mapping_lists[i + 1],
			 data_type[i + 1], f_sa, f_hc, mode);
	pthread_mutex_unlock(&mutex_sp);

	if (data_type[i] == CAL_TYPE) {
	  array_list_clear(mapping_batch->mapping_lists[i], (void *)cal_free);
	} else if (data_type[i] == META_ALIGNMENT_TYPE) {
	  array_list_clear(mapping_batch->mapping_lists[i], (void *)meta_alignment_complete_free);
	} else {
	  array_list_clear(mapping_batch->mapping_lists[i], (void *)alignment_free);
	}
     
	if (data_type[i + 1] == CAL_TYPE) {
	  array_list_clear(mapping_batch->mapping_lists[i + 1], (void *)cal_free);
	} else if (data_type[i + 1] == META_ALIGNMENT_TYPE) {
	  array_list_clear(mapping_batch->mapping_lists[i + 1], (void *)meta_alignment_complete_free);
	} else {
	  array_list_clear(mapping_batch->mapping_lists[i + 1], (void *)alignment_free);
	}
      }

      array_list_free(meta_alignments_list[i], NULL);
      //if (rev_comp[i]) { free(rev_comp[i]); }
      free(scores_ranking[i]);

      array_list_free(meta_alignments_list[i + 1], NULL);
      //if (rev_comp[i + 1]) { free(rev_comp[i + 1]); }
      free(scores_ranking[i + 1]);

      if (post_process_reads[i] == 1 || post_process_reads[i + 1] == 1) {
	//pthread_mutex_lock(&mutex_sp);
	//extern size_t tot_reads_out;
	//tot_reads_out += 2;
	//pthread_mutex_unlock(&mutex_sp);

	fastq_read_free(fq_read0);
	array_list_free(mapping_batch->mapping_lists[i], NULL);

	fastq_read_free(fq_read1);
	array_list_free(mapping_batch->mapping_lists[i + 1], NULL);
      } else {
	//pthread_mutex_lock(&mutex_sp);
	//extern size_t tot_reads_in;
	//tot_reads_in += 2;
	//pthread_mutex_unlock(&mutex_sp);

	array_list_set_flag(1, mapping_batch->mapping_lists[i]);
	mapping_batch->mapping_lists[num_new_reads++] = mapping_batch->mapping_lists[i];
	array_list_insert(fq_read0, new_fq_batch);

	array_list_set_flag(1, mapping_batch->mapping_lists[i + 1]);
	mapping_batch->mapping_lists[num_new_reads++] = mapping_batch->mapping_lists[i + 1];
	array_list_insert(fq_read1, new_fq_batch);
      }
    }
  } else {
    //Single end mode
    for (i = 0; i < num_reads; i++) {
      array_list_free(meta_alignments_list[i], NULL);
      //if (rev_comp[i]) { free(rev_comp[i]); }
      free(scores_ranking[i]);

      fq_read = array_list_get(i, mapping_batch->fq_batch);      
      if (data_type[i] == CAL_TYPE || 
	  data_type[i] == META_ALIGNMENT_TYPE) {
	int mode;
	if (data_type[i] == CAL_TYPE) {
	  mode = 0;
	} else {
	  mode = 1;
	}
	pthread_mutex_lock(&mutex_sp);
	file_write_items(fq_read, mapping_batch->mapping_lists[i],
			 data_type[i], f_sa, f_hc, mode);
	pthread_mutex_unlock(&mutex_sp);
	
	if (data_type[i] == CAL_TYPE) {
	  array_list_clear(mapping_batch->mapping_lists[i], (void *)cal_free);
	} else if (data_type[i] == META_ALIGNMENT_TYPE) {
	  array_list_clear(mapping_batch->mapping_lists[i], (void *)meta_alignment_complete_free);
	}
      }
      //array_list_clear(alignments_list, alignment_free);      
      if (post_process_reads[i] == 0) {
	mapping_batch->mapping_lists[num_new_reads++] =  mapping_batch->mapping_lists[i];
	array_list_insert(fq_read, new_fq_batch);
      } else {
	fastq_read_free(fq_read);
	array_list_free(mapping_batch->mapping_lists[i], NULL);
      }      
    }
  }



  //printf("WK_2ph: =============== REPORT ALIGNMENTS END =====================\n");
  free(post_process_reads);
  array_list_free(mapping_batch->fq_batch, NULL);
  mapping_batch->fq_batch = new_fq_batch;
  free(data_type);
  array_list_free(seeds_list, NULL);
  array_list_free(final_positions, NULL);
  free(new_targets);
  sw_multi_output_free(output);
 
  array_list_free(cals_targets, NULL);

  free(scores_ranking);
  //array_list_free(cals_targets, NULL);
  free(meta_alignments_list);
  //free(rev_comp);

  //============================= Delete!!! ===============================
  //For debug... Validate Reads

  //fprintf(stderr, "APPLY RNA LAST END!\n");
  /*for (int i = 0; i < array_list_size(mapping_batch->fq_batch); i++) {
    fq_read = array_list_get(i, mapping_batch->fq_batch);

    array_list_t *alignments_list = mapping_batch->mapping_lists[i];
    if (!array_list_size(alignments_list)) { continue; }

    //printf("--(%i mappings): Validate %s\n", array_list_size(alignments_list), 
    //	   fq_read->id);
    //--- Extract correct position ---//
    int CHROMOSOME;
    size_t START, END;
	
    char *id = fq_read->id;
    int c = 0, len = strlen(id);
    int pos = 0;
    while (c < len && pos < 4) {
      //printf("while pos [%c]\n", id[c]);
      if (id[c++] == '@') { pos++; }
    }
    
    //printf("Actual pos %c\n", id[c]);
    char value[128];
    int p = 0;
    while (id[c] != '@') {
      value[p++] = id[c++];
    }
    c++;
    value[p] = '\0';
    if (strcmp(value, "X") == 0) { CHROMOSOME = 23; }
    else if (strcmp(value, "Y") == 0) { CHROMOSOME = 24; }
    else if (strcmp(value, "MT") == 0) { CHROMOSOME = 25; }
    else { CHROMOSOME = atoi(value); }

    //printf("Actual pos %c\n", id[c]);
    p = 0;
    while (id[c] != '@') {
      //printf("while pos [%c]\n", id[c]);
      value[p++] = id[c++];
    }
    c++;
    value[p] = '\0';
    START = atol(value);


    p = 0;
    while (id[c] != '@') {
      //printf("while pos [%c]\n", id[c]);
      value[p++] = id[c++];
    }
    c++;
    value[p] = '\0';
    END = atol(value);

    //printf("POSITION: [%i:%i-%i]\n", CHROMOSOME, START, END);
    //---                          ---//
    int map = 0;
    for (int j = 0; j < array_list_size(alignments_list); j++) {
      alignment_t *alignment = array_list_get(j, alignments_list);
      if (alignment->chromosome == CHROMOSOME - 1 && 
	  alignment->position >= START && 
	  alignment->position <= END) {
	map = 1;
	break;
      }
    }

    if (!map) {	  
      for (int t = 0; t < array_list_size(alignments_list); t++) {
	alignment_t *alignment = array_list_get(t, alignments_list);
	alignment_print(alignment);
      }	    

      printf("::ERROR:: %s\n", fq_read->id);
      exit(-1);
    }    
    }*/
  //============================= Delete!!! ===============================

  return LAST_RNA_POST_PAIR_STAGE;

  }

int apply_rna_last_hc(sw_server_input_t* input_p, batch_t *batch) {
  int min_intron_size = input_p->min_intron_size;
  size_t max_intron_size = input_p->max_intron_size;
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  sw_optarg_t *sw_optarg = &input_p->sw_optarg;
  genome_t *genome = input_p->genome_p;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  size_t num_targets = mapping_batch->num_targets;
  metaexons_t *metaexons = input_p->metaexons;
  cal_optarg_t *cal_optarg = input_p->cal_optarg_p;
  avls_list_t *avls_list = input_p->avls_list;
  bwt_optarg_t *bwt_optarg = input_p->bwt_optarg_p;
  bwt_index_t *bwt_index = input_p->bwt_index_p;
  linked_list_t *buffer = input_p->buffer;
  linked_list_t *buffer_hc = input_p->buffer_hc;
  int min_score = input_p->min_score;
  //fprintf(stderr, "APPLY RNA LAST START... %i\n", num_reads);

  array_list_t *cals_list, *fusion_cals, *fusion_cals_aux;
  cal_t *cal, *cal_prev, *cal_next, *first_cal, *last_cal;
  fastq_read_t *fq_read;
  linked_list_iterator_t itr;  
  seed_region_t *s, *s_prev, *s_next;
  cigar_code_t *cigar_code, *cigar_code_prev, *cigar_code_aux;
  cigar_code_t *alig_cigar_code;
  cigar_op_t *cigar_op_start, *cigar_op_end, *cigar_op, *cigar_op_prev, *cigar_op_aux;
  //int *new_targets = (int *)calloc(mapping_batch->num_allocated_targets, sizeof(int));
  array_list_t *merge_cals;
  linked_list_t *linked_list;
  seed_region_t *seed_region;

  //float *cals_score = (float *)calloc(100, sizeof(float));
  float score;
  char reference[2048];
  //char reference_prev[2048];
  //char reference_next[2048];
  //char reference_aux[2048];
  char query[2048];
  // char query_revcomp[2048];
  alignment_t *alignment;
  char q[2048];
  char r[2048];

  //char **rev_comp = (char **)calloc(num_reads, sizeof(char *));
  /*
  fusion_coords_t *extrem_coords[2*40*mapping_batch->num_allocated_targets];
  fusion_coords_t *sp_coords[2*40*mapping_batch->num_allocated_targets];
  cigar_code_t *extrem_cigars[2*40*mapping_batch->num_allocated_targets];
  cigar_code_t *sp_cigars[2*40*mapping_batch->num_allocated_targets];
  */
  char *sequence;
  char *query_ref;
  char *quality_map, *query_map;
  //float scores_ranking[mapping_batch->num_allocated_targets][50];
  float *cals_score;
  char cigar_str[1024];

  cigar_op_t *first_op;
  char *match_seq, *match_qual, *optional_fields, *p;
  int match_start, match_len, optional_fields_length, AS;
  float norm_score;

  int delete_not_cigar;
  int padding_left, padding_right, len_query;  
  int num_sw = 0;
  int num_sw_sp = 0;
  int num_sw_ext = 0;
  size_t num_new_targets = 0;
  int coverage;
  int read_length;
  int cal_pos = 0;
  int number_of_best;
  //int num_extrem_ok = 0;
  //int num_sp_ok = 0;

  //Change to report most alignments
  int seed_err_size = 20;
  int n_alignments = 1;
  int target;
  int c, exact_nt;
  size_t genome_start, genome_end, genome_start2, genome_end2, read_start, read_end;
  int flank = 20;
  int n_delete;
  int sw_pos;
  int seeds_nt; 
  int e1;
  int e2;
  int ref_pos;
  int nt_2;
  
  int start_seeding, end_seeding;
  int lim_start, lim_end;

  //int min_intron_size = 40;

  metaexon_t *first_metaexon, *last_metaexon;
  int l_found, r_found; 
  size_t num_cals;
  register size_t t;
  register int i, j;

  int meta_type;
  meta_alignment_t *meta_alignment;
  sw_item_t *sw_item;  
  sw_depth_t sw_depth;
  sw_depth.depth = 0;
  sw_multi_output_t *output = sw_multi_output_new(MAX_DEPTH);
  //array_list_t **sw_items_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  //array_list_t **alignments_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  array_list_t *meta_alignments_list;
  //array_list_t *seeds_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  //array_list_t *final_positions = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED); 
  int make_seeds = 0;

  //array_list_t *cals_targets;// = array_list_new(50, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *alignments_list;
  //int *post_process_reads = (int *)calloc(num_reads, sizeof(int));
  float scores_ranking[num_reads][100];
  int read_nt;
  int seed_size = cal_optarg->seed_size;

  //pthread_mutex_lock(&mutex_sp);
  //extern size_t tot_reads_in;
  //tot_reads_in += num_reads;
  //pthread_mutex_unlock(&mutex_sp);
  
  //Convert CALs in META_ALIGNMENTS 
  for (int i = 0; i < num_reads; i++) {
    fq_read = array_list_get(i, mapping_batch->fq_batch);
    meta_alignments_list = mapping_batch->mapping_lists[i];
    //printf("WK_3ph (%i): %s\n", array_list_size(meta_alignments_list), fq_read->id );
    
    if (array_list_size(meta_alignments_list) <= 0 || 
	array_list_get_flag(meta_alignments_list) != BITEM_META_ALIGNMENTS) { continue; }   
    
    //printf("\tProcess\n");
    for (int m = 0; m < array_list_size(meta_alignments_list); m++) {
      meta_alignment = array_list_get(m, meta_alignments_list);

      fusion_cals = meta_alignment->cals_list;

      cal_t *first_cal = array_list_get(0, fusion_cals);
      cal_t *last_cal = array_list_get(array_list_size(fusion_cals) - 1, fusion_cals);

      //cal_print(first_cal);
      //cal_print(last_cal);

      if (first_cal->strand == 1) {
	//if (!rev_comp[i]) {
	//rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	//strcpy(rev_comp[i], fq_read->sequence);
	//seq_reverse_complementary(rev_comp[i], fq_read->length);
	//}
	query_map = fq_read->revcomp;
	//query_map = rev_comp[i];
      } else {
	query_map = fq_read->sequence;
      }

      s_prev = linked_list_get_first(first_cal->sr_list);
      s_next = linked_list_get_last(last_cal->sr_list);

      if (s_prev->read_start != 0 &&
	  meta_alignment->cigar_left == NULL) {
	/*metaexon_search(first_cal->strand, first_cal->chromosome_id - 1,
	  first_cal->start, first_cal->end, &first_metaexon,
	  metaexons);
	*/
	cigar_code_t *cc_left = fill_extrem_gap(query_map, 
						first_cal,
						FILL_GAP_LEFT,
						genome,
						metaexons, 
						avls_list); 
	//printf("LEFT CIGAR: %s\n", new_cigar_code_string(cc_left));
	//assert(meta_alignment != NULL);
	meta_alignment_insert_cigar(cc_left, CIGAR_ANCHOR_RIGHT, 0, meta_alignment);
      }
	
      if (s_next->read_end != fq_read->length - 1 &&
	  meta_alignment->cigar_right == NULL) {
	/*metaexon_search(last_cal->strand, last_cal->chromosome_id - 1,
	  last_cal->start, last_cal->end, &first_metaexon,
	  metaexons);	    */
	cigar_code_t *cc_right = fill_extrem_gap(query_map, 
						 last_cal,
						 FILL_GAP_RIGHT,
						 genome,
						 metaexons, 
						 avls_list); 	  
	//printf("RIGHT CIGAR: %s\n", new_cigar_code_string(cc_right));
	meta_alignment_insert_cigar(cc_right, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
      }

      //SW for complete meta-alignments extrems
      //meta_alignments_list = mapping_batch->mapping_lists[i];
      //for (int j = 0; j < array_list_size(meta_alignments_list); j++) {
      //meta_alignment = array_list_get(j, meta_alignments_list);
      //cal = array_list_get(0, meta_alignment->cals_list);
      //cal_t *first_cal = array_list_get(0, fusion_cals);
      //cal_t *last_cal = array_list_get(array_list_size(fusion_cals) - 1, fusion_cals);

      /*array_list_t *fusion_cals = meta_alignment->cals_list;
	printf("==== FUSION CALS PROCESS ====\n");
	for (int t = 0; t < array_list_size(fusion_cals); t++) {
	cal_t *cal_aux = array_list_get(t, fusion_cals);
	cal_print(cal_aux);
	}	
      */    
      //printf("Meta - %i:%lu-%lu\n", cal->chromosome_id, cal->start, cal->end);
      //if (cal->strand == 1) {
      //if (!rev_comp[i]) {
      //rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
      //strcpy(rev_comp[i], fq_read->sequence);
      //seq_reverse_complementary(rev_comp[i], fq_read->length);
      //}
      //query_map = rev_comp[i];
      //} else {
      //query_map = fq_read->sequence;
      //}
      //s_prev = linked_list_get_first(first_cal->sr_list);
      //s_next = linked_list_get_last(last_cal->sr_list);

      const int flank_s = 5;
      
      if (meta_alignment->cigar_left == NULL && 
	  s_prev->read_start  != 0) {
	//seed_region_t *seed_reg = linked_list_get_first(cal->sr_list);
	if (s_prev->read_start >= flank_s) {
	  //seed_reg->read_start += flank_s;
	  //seed_reg->genome_start += flank_s;
	  //cal->start += flank_s;	  
	  //cigar_code_delete_nt(flank_s, 0, cal->info);	  
	  //SW
	  genome_start = s_prev->genome_start - s_prev->read_start ;
	  genome_end   = s_prev->genome_start - 1;
	  genome_read_sequence_by_chr_index(r, 0, 
					    first_cal->chromosome_id - 1,
					    &genome_start, &genome_end,
					    genome);    
	  memcpy(q, query_map, s_prev->read_start);
	  q[s_prev->read_start] = '\0';	  
	  //printf("query     L ::: %s\n", q);
	  //printf("reference L ::: %s\n", r);
	  //New sw item. Storing data...
	  sw_item = sw_item_new(EXTREM_SW_LEFT, i, j, j,
				first_cal, first_cal, 
				meta_alignment, NULL, 
				NULL, NULL);	      
	  //Insert item... and process if depth is full
	  sw_depth_insert(q, r, sw_item,
			  sw_optarg, output,
			  avls_list, metaexons, 
			  &sw_depth, genome, min_intron_size); 
	} else {
	  cigar_code = cigar_code_new();
	  cigar_code_append_new_op(s_prev->read_start, 'M', cigar_code);
	  meta_alignment->cigar_left = cigar_code;
	}
      }

      if (meta_alignment->cigar_right == NULL &&
	  s_next->read_end != fq_read->length - 1) {
	//seed_region_t *seed_reg = linked_list_get_last(s_next->sr_list);
	if ((fq_read->length - 1) - s_next->read_end >= flank_s) {
	  //cal = array_list_get(array_list_size(meta_alignment->cals_list) - 1, meta_alignment->cals_list);
	  //seed_reg->read_end -= flank_s;
	  //seed_reg->genome_end -= flank_s;
	  //cal->end -= flank_s;

	  //cigar_code_delete_nt(flank_s, 1, cal->info);
	  
	  //if (seed_reg->read_end <= fq_read->length - 5) {
	  //SW
	  int r_gap = fq_read->length - s_next->read_end - 1;
	  
	  //printf("r_gap = %i, genome_end = %lu, read_end = %i\n", 
	  //	 r_gap, seed_reg->genome_end + 1, seed_reg->read_end + 1);

	  //cal_print(last_cal);
	  
	  genome_start = s_next->genome_end + 1;
	  genome_end   = genome_start + r_gap - 1;	  
	  genome_read_sequence_by_chr_index(r, 0,
					    last_cal->chromosome_id - 1,
					    &genome_start, &genome_end,
					    genome);    
	  memcpy(q, query_map + s_next->read_end + 1, r_gap);
	  q[r_gap] = '\0';
	  //printf("query     R ::: %s\n", q);
	  //printf("reference R (%lu-%lu)::: %s\n", genome_start, genome_end, r);
	  
	  //New sw item. Storing data...
	  sw_item = sw_item_new(EXTREM_SW_RIGHT, i, j, j,
				last_cal, last_cal, 
				meta_alignment, NULL, 
				NULL, NULL);	      
	  //Insert item... and process if depth is full
	  sw_depth_insert(q, r, sw_item,
			  sw_optarg, output,
			  avls_list, metaexons, 
			  &sw_depth, genome, min_intron_size); 
	} else {
	  cigar_code = cigar_code_new();
	  cigar_code_append_new_op(fq_read->length - s_next->read_end - 1, 'M', cigar_code);
	  meta_alignment->cigar_right = cigar_code;
	}
      }
      //fprintf(stderr, "WK_3ph: %s FINISH\n", fq_read->id);
    }
  }

  sw_depth_process(sw_optarg, output, 
		   &sw_depth, avls_list,
		   metaexons, SW_FINAL,
		   genome, min_intron_size);

  const int MIN_SCORE_CAL = 50;
  //printf("=============== REPORT ALIGNMENTS =====================\n");
  //fprintf(stderr, "REPORT META ALIGNMENTS\n");
  //fprintf(stderr, "WK_3ph: =============== REPORT ALIGNMENTS (%i) =====================\n", num_reads);
  for (int i = 0; i < num_reads; i++) {
    fq_read = array_list_get(i, mapping_batch->fq_batch);
    //printf("WK_3ph-Meta_Report (%i):  %s >>>>\n", array_list_size(mapping_batch->mapping_lists[i]),
    //	   fq_read->id);

    meta_alignments_list = mapping_batch->mapping_lists[i];    

    if (array_list_size(meta_alignments_list) <= 0 || 
	array_list_get_flag(meta_alignments_list) != BITEM_META_ALIGNMENTS) { 
      array_list_set_flag(1, mapping_batch->mapping_lists[i]);
      continue;
    }   

    assert(fq_read->id != NULL);
    //fprintf(stderr, ".( %i )Read %i: %s .\n", array_list_size(meta_alignments_list[i]), i, fq_read->id);
    //printf(".( %i )Read %i: %s .\n", array_list_size(meta_alignments_list[i]), i, fq_read->id);
    
    char query[2048];
    char quality[2048];
    int map = 0;

    for (int m = array_list_size(meta_alignments_list) - 1; m >= 0; m--) {
      meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list);
      meta_alignment_close(meta_alignment);

      if (meta_alignment_get_status(meta_alignment) != META_CLOSE) { 
	array_list_remove_at(m, meta_alignments_list);
	//printf("Meta alignment Closed? [-NOT CLOSED-]\n");
	meta_alignment_complete_free(meta_alignment);
	continue;
      } else {
	//printf("Meta alignment Closed? [-CLOSED-]\n");
      }
      meta_alignment_calculate_score(meta_alignment);
    }

    if (array_list_size(meta_alignments_list) <= 0) { continue; }

    alignments_list = array_list_new(array_list_size(mapping_batch->mapping_lists[i]), 
				     1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    
    size_t start_mapping;
    int n_report = array_list_size(meta_alignments_list);//array_list_size(meta_alignments_list) >= 2 ? 2 : array_list_size(meta_alignments_list);
    //int n_report = array_list_size(meta_alignments_list) >= 5 ? 5 : array_list_size(meta_alignments_list);

    for (int m = 0; m < n_report; m++) {
      meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list);
      map = 1;
      optional_fields_length = 0;
      optional_fields = NULL;
      
      first_cal = meta_alignment_get_first_cal(meta_alignment);
      if (first_cal->strand == 1) {
	//if (!rev_comp[i]) {
	//rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	//strcpy(rev_comp[i], fq_read->sequence);
	//seq_reverse_complementary(rev_comp[i], fq_read->length);
	//}
	query_map = fq_read->revcomp;
	//query_map = rev_comp[i];
      } else {
	query_map = fq_read->sequence;
      }

      cigar_code = meta_alignment->cigar_code;	
      assert(cigar_code != NULL);       	    

      int dsp = 0;
      start_mapping = first_cal->start;
      s_prev = linked_list_get_first(first_cal->sr_list);
      if (meta_alignment->cigar_left == NULL && 
	  s_prev->read_start > 0) {
	cigar_op_t *op = cigar_op_new(s_prev->read_start, 'H');
	array_list_insert_at(0, op, cigar_code->ops);
      } else {
	if (meta_alignment->cigar_left != NULL ) {
	  //printf("Cigar %s\n", new_cigar_code_string(meta_alignment->cigar_left));
	  cigar_code_t *cigar_code = meta_alignment->cigar_left;
	  for (int c = 0; c < cigar_code->ops->size; c++) {
	    cigar_op_t *op = array_list_get(c, cigar_code->ops);
	    if (op->name == 'M' ||
		op->name == 'N' ||
		op->name == 'D') {
	      dsp += op->number;
	    }
	  }
	}
      }
	
      start_mapping -= dsp;
      //printf(" new_start = %lu\n", start_mapping);

      last_cal = array_list_get(array_list_size(meta_alignment->cals_list) - 1, meta_alignment->cals_list);
      s_next = linked_list_get_last(last_cal->sr_list);
      //printf("LAST SEED %i, report %i, %p, %p\n", s_next->read_end, fq_read->length - s_next->read_end - 1, s_prev, s_next);

      if (meta_alignment->cigar_right == NULL && 
	  s_next->read_end < fq_read->length - 1) {
	//printf("Cigar NULL R\n");
	cigar_op_t *op = cigar_op_new(fq_read->length - s_next->read_end - 1, 'H');
	array_list_insert(op, cigar_code->ops);
      }
      

      //Cigar  complete! filter bad alignments by score!
      //-------------------------------------------------

      //printf("%s: CIGAR : %s\n", fq_read->id, new_cigar_code_string(cigar_code));

      int n_m = 0, n_d = 0, n_i = 0;
      int tot_dist = cigar_code->distance;

      cigar_op_t *op_a = array_list_get(0, cigar_code->ops);
      if (op_a->name == 'H') { tot_dist += op_a->number; }

      op_a = array_list_get(array_list_size(cigar_code->ops) - 1, cigar_code->ops);
      if (op_a->name == 'H') { tot_dist += op_a->number; }

      for (int c = 0; c < cigar_code->ops->size; c++) {
	cigar_op_t *op = array_list_get(c, cigar_code->ops);
	if (op->name == 'M') { n_m += op->number; }
	else if (op->name == 'D') { n_d += op->number; }
	else if (op->name == 'I') { n_i += op->number; }
      }

      if (tot_dist > fq_read->length / 2) { continue; }

      //printf("(num ops %i) (dist tot %i | distance %i) NUM MATCHES = %i, NUM_DEL = %i, NUM_INSERT = %i\n", 
      //     cigar_code->ops->size, tot_dist, cigar_code->distance, n_m, n_d, n_i);

      //-------------------------------------------------


      if (first_cal->strand == 1) {
	//query_map = rev_comp[i];
	query_map = fq_read->revcomp;
      } else {
	query_map = fq_read->sequence;
      }
		
      int h_left, h_right;
      cigar_op_t *first_op = array_list_get(0, cigar_code->ops);
      if (first_op->name == 'H') {	  
	h_left = first_op->number;
      } else  {
	h_left = 0;
      }
	
      cigar_op_t *last_op = array_list_get(array_list_size(cigar_code->ops) - 1, cigar_code->ops);
      if (last_op->name == 'H') { 
	h_right = last_op->number;
      } else  {
	h_right = 0;
      }

      if (!cigar_code_validate_(fq_read, cigar_code)) {
	meta_alignment_complete_free(meta_alignment);
	continue;
      }
      
      int read_score = cigar_code_score(cigar_code, fq_read->length);
      
      if (read_score < min_score) { 
	meta_alignment_complete_free(meta_alignment);
	continue;
      }

      //pthread_mutex_unlock(&mutex_sp);

      //printf("FINAL H_LEFT = %i, H_RIGHT = %i\n", h_left, h_right);
      //printf("FINAL H_LEFT = %i, H_RIGHT = %i\n", h_left, h_right);

	
      //if (h_left > fq_read->length || h_left < 0) { exit(-1); }
      //if (h_left + h_right >= fq_read->length) { continue; }
	
      int len_read = fq_read->length - h_left - h_right;
      //printf("(REAL %i): %s (%i)\n", len_read, query, s);
	
      //printf("* * * %s * * *\n", fq_read->id);
      //fprintf(stderr, "* * * (%i) M E T A    A L I G N M E N T    R E P O R T    %s* * *\n", strlen(query), new_cigar_code_string(cigar_code));
      alignment = alignment_new();
	
      //int  header_len = strlen(fq_read->id); 
      //char header_id[header_len + 1];
      //get_to_first_blank(fq_read->id, header_len, header_id);
      //char *header_match = (char *)malloc(sizeof(char)*header_len);
      //if (header_match == NULL) { exit(-1); }
      //memcpy(header_match, header_id, header_len);
		
      memcpy(query, &query_map[h_left], len_read);
      query[len_read] = '\0';
	

      memcpy(quality, &fq_read->quality[h_left], len_read);
      quality[len_read] = '\0';
      
      /*pthread_mutex_lock(&mutex_sp);
        float score_tmp = len_read * 0.5 - cigar_code->distance * 0.4;
        float score_map = (score_tmp * 254) / (len_read * 0.5);
	extern float min_score;
        if (score_map < min_score) {
          min_score = score_map;
          printf("Read (score_tmp = %0.2f , score_map = %0.2f, distance = %i, len_read = %i, %0.2f): %s\n",
                 score_tmp, score_map, min_score, cigar_code->distance, len_read, fq_read->id);
        }
        extern float tot_score;
	extern size_t tot_rep;
        tot_score += score_map;
        tot_rep++;
        pthread_mutex_unlock(&mutex_sp);
      */

      alignment_init_single_end(strdup(fq_read->id),
				strdup(query)/*match_seq*/,
				strdup(quality)/*match_qual*/,
				first_cal->strand, first_cal->chromosome_id - 1, start_mapping - 1,
				strdup(new_cigar_code_string(cigar_code))/*strdup(cigar_fake)*/,
				cigar_code_get_num_ops(cigar_code)/*1*/,
				cigar_code_score(cigar_code, fq_read->length),
				1, (array_list_size(meta_alignments_list) < 1),
				cigar_code->distance, NULL, alignment);
      
      //printf("Report CIGAR OK!\n");	
      array_list_insert(alignment, alignments_list);	
      meta_alignment_complete_free(meta_alignment);
      
    }
      //}

    array_list_free(mapping_batch->mapping_lists[i], NULL);
    mapping_batch->mapping_lists[i] = alignments_list;
    array_list_set_flag(1, mapping_batch->mapping_lists[i]);

    size_t n_alignments = array_list_size(alignments_list);
    int final_distance;
    for (size_t a = 0; a < n_alignments; a++) {
      alignment_t *alignment = (alignment_t *)array_list_get(a, alignments_list);
      
      // set optional fields                                                             	
      optional_fields_length = 100;
      optional_fields = (char *) calloc(optional_fields_length, sizeof(char));
      
      p = optional_fields;
      AS = alignment->map_quality;
      alignment->map_quality = 255;
      //AS = (int) 254;
      
      sprintf(p, "ASi");
      p += 3;
      memcpy(p, &AS, sizeof(int));
      p += sizeof(int);
      
      sprintf(p, "NHi");
      p += 3;
      memcpy(p, &n_alignments, sizeof(int));
      p += sizeof(int);
      
      sprintf(p, "NMi");
      p += 3;
      final_distance = alignment->optional_fields_length;
      memcpy(p, &final_distance, sizeof(int));
      p += sizeof(int);	

      alignment->optional_fields_length = optional_fields_length;
      alignment->optional_fields = optional_fields;
    }
    //if (rev_comp[i]) { free(rev_comp[i]); }
  }

  sw_multi_output_free(output);
  //free(rev_comp);
  
  return LAST_RNA_POST_PAIR_STAGE;
  
}
