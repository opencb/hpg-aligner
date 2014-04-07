/*
 * rna_server.c
 *
 *  Created on: Jan 20, 2012
 *  Last modified: Dec 20, 2013
 *  Author: Hector Martinez
 */

#ifndef RNA_SERVER_H
#define RNA_SERVER_H

#include "breakpoint.h"
#include "sw_server.h"
#include "buffers.h"
#include "sw.h"
#include "cal_seeker.h"
#include  "preprocess_rna.h"

#include "containers/list.h"
#include "commons/commons.h"
#include "commons/system_utils.h"
#include "bioformats/bam/alignment.h"
#include "aligners/bwt/bwt.h"
#include "aligners/bwt/genome.h"
#include "aligners/sw/macros.h"

#define FILL_GAP_LEFT  0
#define FILL_GAP_RIGHT 1

#define META_ALIGNMENT_NONE   0
#define META_ALIGNMENT_MIDDLE 1
#define META_ALIGNMENT_LEFT   2
#define META_ALIGNMENT_RIGHT  3

//#include "list.h"
//#include "genome.h"

//===============================================================================================================

typedef struct info_sp {
  size_t l_genome_start;
  size_t l_genome_end;
  size_t r_genome_start;
  size_t r_genome_end;
} info_sp_t;

info_sp_t* sw_reference_splice_junction(cal_t *cal_prev, cal_t *cal_next,
					char *query_map, genome_t *genome,
					char *q, char *r);

/**
 * @brief Structure for store splice junctions information.
 * 
 * Structure for store the strand, extend and exact coordiantes for each splice junction found.  
 */
typedef struct sp_data{
  size_t start_sp; /**< exact start genome coordinate */
  size_t end_sp;  /**< exact end genome coordinate */
  short int strand_sp;  /**< genome strand */
  size_t start_extend_sp; /**< extend start genome coordinate */ 
  size_t end_extend_sp;  /**< extend end genome coordinate */
}sp_data_t;

/**
 * @brief Structure for store CAL data information.
 * 
 * Structure for store genome coordinates and string coordinates of CALs fusion 
 */
typedef struct cal_fusion_data{
  unsigned int  id;   /**< id of CAL */
  size_t  genome_start; /**< genome start position */
  size_t  genome_end;   /**< genome end position */
  unsigned char genome_strand;   /**< genome strand position */
  unsigned int  genome_chromosome;  /**< genome chromosome */
  unsigned int  fusion_start;  /**< fusion string start position */
  unsigned int  fusion_end;   /**< fusion string end position */
}cal_fusion_data_t;


/**
 * @brief  Initializer for the @a cal_fusion_data_t structure.
 * @param  id of CAL
 * @param  start genome start position
 * @param  end genome end position
 * @param  strand genome strand position
 * @param  chromosome genome chromosome positio
 * @param  fusion_start fusion string start position
 * @param  fusion_end fusion string end position

 *  Initialize all @a cal_fusion_data_t fields with the information CALs.
 */
void cal_fusion_data_init(unsigned int id, size_t start, size_t end, unsigned int strand, unsigned int chromosome, unsigned int fusion_start, unsigned int fusion_end, cal_fusion_data_t *cal_fusion_data_p);


/**
 * @brief Structure for allocate all @a cal_fusion_data_t structures.
 * 
 * Structure for store the CAL data needed for process CALs fusion after smith-waterman process ends.  
 */
typedef struct allocate_fusion_data{
  cal_fusion_data_t *allocate_data; /**< array to store all CAL data structures */ 
  unsigned int cal_number; /**< number of CAL data structures are allocated in the array */
}allocate_fusion_data_t;

/**
 * @brief  Process Smith-Waterman algorithm for all reads with their CALs and then search splice junctions.
 * @param  input_p input_p all configuration values and some data structures needed by the server
 * @param  chromosome_avls_p array that allocates the avls trees on stored all splice junctions
 *
 * For each read process their CALs with Smith-Waterman algorithm and then search for each output splice junctions.
 * For it, search big gaps in Smith-Waterman output sequence and then search cannonical splice junctions (GT-AG and CT-AC).
 * Finally generates cigar string with the mappings results and stores all splice junctions in avl trees.
 */
void rna_server_omp_smith_waterman(sw_server_input_t* input_p, allocate_splice_elements_t *chromosome_avls_p);

/**
  * @brief Structure for store each cigar operation.
  *
  * Structure for store the type of cigar operation and the value of this.
  */
typedef struct cigar_data{
  unsigned char type;  /**< type of cigar operation */
  int value;  /**< value of cigar operation */
}cigar_data_t;

/**
 * @brief  Return the description of one cigar operation.
 * @param  status code of one cigar operation
 * @return return the description of cigar operation
 * 
 * Return the description of the cigar operation input.
 */
//char* cigar_automata_status(unsigned char status);

/**
 * @brief  Insert new cigar operation to @a cigar_data_t.
 * @param  cigar_p array for store all cigar operations
 * @param  max_len maximum number for allocate cigar operations
 * @param  pos actual array position 
 * @param  status type of new cigar operation
 * @param  number value of new cigar operation 
 * 
 * Insert one new cigar operation to the allocate array.
 */
//void cigar_generator(cigar_data_t *cigar_p, unsigned int *max_len, unsigned int *pos, unsigned int status, unsigned int *number);

/**
 * @brief  Process Smith-Waterman outputs for search splice junctions and generate cigar string.
 * @param  sw_output_p Smith-Waterman output
 * @param  sw_input_p Smith-Waterman input
 * @param  depth number of Smith-Waterman processed
 * @param  depth_cal_fusion_p data allocate for each CAL fusioned
 * @param  chromosome_avls_p avls trees for store all splice junctions found
 * @param  sw_channels_p data stored for each level of depth for Smith-Waterman process
 * @param  sw_batch_p structure that stores data for each read
 * @param  write_list_p list for store write batches
 * @param  write_batch_p batch of data that will be written in a file
 * @param  write_size size of write batch
 * @param  sw_id Smith-Waterman process id
 * @param  sw_no_valids number of Smith-Watermans no valids
 * @param  min_score minimum score for process Smith-Waterman output
 * @param  genome_p structure that allocates all genome
 * @param  min_intron_length minimum intron length for store it
 *
 * 
 * After observe if the Smith-Waterman score output is major that @a min_score 
 * process Smith-Waterman output for search big gaps, search splice junctions and 
 * generate cigar string. 
 */
void search_splice_junctions_sw_output(sw_simd_input_t* sw_input_p, 
				       sw_simd_output_t* sw_output_p, 
				       unsigned int depth, 
				       allocate_fusion_data_t *depth_cal_fusion_p,
				       allocate_splice_elements_t *chromosome_avls_p, 
				       sw_channel_t *sw_channels_p, 
				       mapping_batch_t *mapping_batch_p, 
				       unsigned int sw_id, 
				       size_t *sw_no_valids, 
				       float min_score, 
				       genome_t *genome_p, 
				       int min_intron_length, 
				       float gap_open, 
				       float gap_extend,
				       float match);

int apply_sw_rna(sw_server_input_t* input_p, batch_t *batch);

int apply_rna_last(sw_server_input_t* input_p, batch_t *batch);

int apply_rna_last_hc(sw_server_input_t* input_p, batch_t *batch);

int search_simple_splice_junction(seed_region_t *s_prev, seed_region_t *s_next,
				  int chromosome_id, int strand, 
				  char *sequence, genome_t *genome, 
				  size_t *sp_start, size_t *sp_end,
				  int *sp_type,
				  int *distance);

int search_simple_splice_junction_semi_cannonical(seed_region_t *s_prev, seed_region_t *s_next,
						  int chromosome_id, int strand, 
						  char *sequence, genome_t *genome, 
						  size_t *sp_start, size_t *sp_end,
						  int *sp_type,
						  int *distance);

cigar_code_t *search_double_anchors_cal(char *query_map,
					cal_t *first_cal, 
					cal_t *last_cal,
					metaexons_t *metaexons, 
					genome_t *genome,
					fastq_read_t *fq_read,
					int *type, 
					avls_list_t *avls_list);

cigar_code_t *search_right_single_anchor(int gap_close, 
					 cal_t *cal,
					 int filter_pos, 
					 array_list_t *left_breaks,
					 char *query_map, metaexons_t *metaexons,
					 genome_t *genome,
					 avls_list_t *avls_list);

cigar_code_t *search_left_single_anchor(int gap_close, 
					cal_t *cal,
					int filter_pos, 
					array_list_t *right_breaks,
					char *query_map,
					metaexons_t *metaexons,
					genome_t *genome,
					avls_list_t *avls_list);

cigar_code_t *fill_extrem_gap(char *query, 
			      cal_t *cal,
			      int type,
			      genome_t *genome,
			      metaexons_t *metaexons, 
			      avls_list_t *avls_list);


cigar_code_t *meta_alignment_fill_extrem_gap(char *query, 
					     cal_t *cal,
					     int type,
					     genome_t *genome,
					     metaexons_t *metaexons, 
					     avls_list_t *avls_list);

#endif
