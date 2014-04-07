#ifndef BUFFERS_H
#define BUFFERS_H

#include "containers/array_list.h"
#include "containers/list.h"
#include "aligners/bwt/bwt.h"
#include "bioformats/bam/alignment.h"
#include "aligners/bwt/genome.h"
#include "timing.h"
#include "statistics.h"
#include "commons/log.h"

#include "breakpoint.h"

#include "bs/array_list_bs.h"

//#include "bwt_server.h"
//#include "rna/rna_server.h"

//-----------------------------------------------

#define ALIGNMENT_TYPE       0
#define CAL_TYPE             1
#define META_ALIGNMENT_TYPE  2

//====================================================================================
//  Buffer RNA Vars
//====================================================================================

#define BITEM_NO_CALS          0
#define BITEM_SINGLE_ANCHORS   1
#define BITEM_CALS             2
#define BITEM_META_ALIGNMENTS  3

//====================================================================================


//====================================================================================
//  Workflow Vars
//====================================================================================

//================================= NORMAL WORKFLOW ==================================

//-------- DEFINE WORKFLOW COMMON VARS -----------

#define BWT_STAGE               0
#define CONSUMER_STAGE         -1


// added by PP
#define BS_BWT_STAGE            0
//#define BS_SEEDING_STAGE        1
#define BS_CAL_STAGE            1
#define BS_PRE_PAIR_STAGE       2
#define BS_SW_STAGE             3
#define BS_POST_PAIR_STAGE      4
#define BS_STATUS_STAGE         5

//--------  DEFINE WORKFLOW DNA VARS  -----------


#define CAL_STAGE               1
#define PRE_PAIR_STAGE          2
#define SW_STAGE                3
#define DNA_POST_PAIR_STAGE     4


//--------  DEFINE WORKFLOW RNA VARS   -----------

#define RNA_CAL_STAGE           1
#define RNA_STAGE               2
#define RNA_POST_PAIR_STAGE     3

//================================= LAST WORKFLOW ==================================

#define LAST_RNA_POST_PAIR_STAGE     1

//====================================================================================
//  globals
//====================================================================================

#define NOT_ANCHORS         0
#define SINGLE_ANCHORS      1
#define DOUBLE_ANCHORS      2
#define ALIGNMENTS_FOUND    3
#define ALIGNMENTS_EXCEEDED 4

//------------------------------------------------------------------------------------

#define DNA_MODE           0
#define RNA_MODE           1
// added by PP
#define BS_MODE            2

//------------------------------------------------------------------------------------

#define DNA_FLAG           1
#define RNA_FLAG           2
#define SINGLE_END_FLAG    4
#define PAIRED_END_FLAG    8
#define MATE_PAIR_FLAG    16
#define PAIR1_FLAG        32
#define PAIR2_FLAG        64
#define WRITE_ITEM_FLAG  128
#define SW_ITEM_FLAG     256
// added by PP
#define BS_FLAG          512

//------------------------------------------------------------------------------------

#define SINGLE_END_MODE 0
#define PAIRED_END_MODE 1
#define MATE_PAIR_MODE  2

//------------------------------------------------------------------------------------

#define PAIR_1   1
#define PAIR_2   2
#define PAIR_1_2 3

//------------------------------------------------------------------------------------

#define UNKNOWN_ITEM    0
#define READ_ITEM       1
#define KL_ITEM         2
#define SEED_ITEM       3
#define SPLIT_KL_ITEM   4
#define SW_ITEM         5
#define WRITE_ITEM      6

//------------------------------------------------------------------------------------

#define MISMATCH_FLAG        0
#define MATCH_FLAG           1
#define SPLICE_EXACT_FLAG    2
#define SPLICE_EXTEND_FLAG   3

//------------------------------------------------------------------------------------

#define NORMAL_MODE 0
#define SEED_MODE   1

//------------------------------------------------------------------------------------

#define FASTQ_READER         	0
#define BWT_SERVER	   	1
#define REGION_SEEKER		2
#define CAL_SEEKER	   	3
#define RNA_PREPROCESS	    	4
#define RNA_SERVER	    	5
#define BAM_WRITER         	6
#define TOTAL_TIME         	7

//------------------------------------------------------------------------------------

#define MAX_READ_MAPPINGS 100 

//------------------------------------------------------------------------------------

#define UNKNWON_ACTION 0
#define BWT_ACTION     1
#define SEEDING_ACTION 2
#define CAL_ACTION     3
#define PAIR_ACTION    4
#define SW_ACTION      5

#define NUM_STRANDS 2

//====================================================================================
//  SPLICE JUNCTION TYPE
//====================================================================================

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


//====================================================================================
//  structures and prototypes
//====================================================================================

bam_header_t *create_bam_header_by_genome(genome_t *genome);


/**
 * @brief Structure for store in files all the process data.
 * 
 * Structure for store in files all data process by each pipeline phase.
 */
typedef struct write_batch {
  unsigned char flag;           /**< type of data stored*/
  unsigned int size;            /**< actual size of the batch (in bytes)*/
  unsigned int allocated_size;  /**< maximum size of the batch (in bytes)*/
  void* buffer_p;               /**< buffer to store data*/
} write_batch_t;

//-------------------------------------------------------------------------------------

/**
 * @brief  Constructor for the @a write_batch_t structure.
 * @param  allocate_size maximum size of the batch (in bytes)
 * @param  flag type of data stored
 * @return Pointer to the new structure.
 * 
 * Constructor function that allocates memory for
 * the batch_writer structure and initialize it.
 */
write_batch_t* write_batch_new(unsigned int allocate_size, unsigned char flag);

/**
 * @brief  Destrcutor for the @a write_batch_t structure.
 * @param  write_batch_p[out] pointer to the structure to free
 * 
 * @a write_batch_t destructor that frees the memory previosly
 * allocated by constructor @a write_batch_new
 */
void write_batch_free(write_batch_t* write_batch_p);

//====================================================================================

typedef struct region_batch {
  array_list_t **allocate_mapping_p;
  fastq_batch_t *unmapped_batch_p;
} region_batch_t;

void region_batch_init(array_list_t **allocate_mapping_p, fastq_batch_t *unmapped_batch_p, 
		       region_batch_t *region_batch_p);
void region_batch_free(region_batch_t *allocate_regions_p);

//====================================================================================
/*
typedef struct sw_batch {
  unsigned int num_reads;
  array_list_t **allocate_cals_p;
  fastq_read_t **allocate_reads_p;
} sw_batch_t;

void sw_batch_init(unsigned int num_reads, array_list_t **allocate_cals_p, 
		   fastq_read_t **allocate_reads_p, sw_batch_t *sw_batch_p);
void sw_batch_free(sw_batch_t *sw_data_p);
*/
//====================================================================================

typedef struct report_optarg {
  int all;
  int n_best;
  int n_hits;
  int only_paired;
  int best;
} report_optarg_t;

report_optarg_t *report_optarg_new(int all, int n_best, int n_hits, int only_paired, int best);

void report_optarg_free(report_optarg_t *p);

//====================================================================================

typedef struct pair_mng {
  int pair_mode;
  size_t min_distance;
  size_t max_distance;
  int report_only_paired;
} pair_mng_t;

pair_mng_t *pair_mng_new(int pair_mode, size_t min_distance, 
			 size_t max_distance, int report_only_upaired);

void pair_mng_free(pair_mng_t *p);

//=====================================================================================

// Added by PP
//====================================================================================

typedef struct bs_context {
  size_t CpG_methyl;                 /**< Partial Counter for methylated Cytosines in CpG context   */
  size_t CpG_unmethyl;               /**< Partial Counter for unmethylated Cytosines in CpG context */
  size_t CHG_methyl;                 /**< Partial Counter for methylated Cytosines in CHG context   */
  size_t CHG_unmethyl;               /**< Partial Counter for unmethylated Cytosines in CpG context */
  size_t CHH_methyl;                 /**< Partial Counter for methylated Cytosines in CHH context   */
  size_t CHH_unmethyl;               /**< Partial Counter for unmethylated Cytosines in CpG context */
  size_t MUT_methyl;                 /**< Partial Counter for mutated Cytosines                     */
  size_t num_bases;                  /**< Partial Counter for number of bases in the batch          */
  array_list_t *context_CpG;             /**< Array with the sequences from CpG context to write */
  array_list_t *context_CHG;             /**< Array with the sequences from CHG context to write */
  array_list_t *context_CHH;             /**< Array with the sequences from CHH context to write */
  array_list_t *context_MUT;             /**< Array with the sequences from mutations to write   */

  array_list_bs_t *context_bs_CpG;             /**< Array with the sequences from CpG context to write */
  array_list_bs_t *context_bs_CHG;             /**< Array with the sequences from CHG context to write */
  array_list_bs_t *context_bs_CHH;             /**< Array with the sequences from CHH context to write */
  array_list_bs_t *context_bs_MUT;             /**< Array with the sequences from mutations to write   */
} bs_context_t;

bs_context_t *bs_context_new(size_t num_reads);
void bs_context_free(bs_context_t *bs_context);

void bs_context_init(bs_context_t * bs_context, size_t num_reads);

//====================================================================================

typedef struct mapping_batch {
  int action;
  size_t num_targets;
  size_t num_extra_targets;
  size_t num_allocated_targets;
  size_t num_to_do;
  unsigned char extra_stage_do;
  unsigned char was_process;

  size_t num_gaps;
  size_t num_sws;
  size_t num_ext_sws;

  unsigned char *extra_stage_id;
  array_list_t *fq_batch;
  size_t *targets;
  size_t *extra_targets;
  array_list_t **mapping_lists;
  char *status; 
  pair_mng_t *pair_mng;
  array_list_t **old_mapping_lists;
  unsigned char *bwt_mappings;

  size_t *histogram_sw;

  // bs handling
  size_t num_targets2;
  size_t num_to_do2;
  size_t *targets2;

  array_list_t **mapping_lists2;
  array_list_t *CT_fq_batch;
  array_list_t *GA_fq_batch;

  array_list_t *CT_rev_fq_batch;
  array_list_t *GA_rev_fq_batch;

  array_list_t *bs_status;
  bs_context_t *bs_context;
  //  bs_context_t bs_context;
} mapping_batch_t;

mapping_batch_t *mapping_batch_new_2(size_t num_reads, array_list_t *fq_batch, pair_mng_t *pair_mng);
mapping_batch_t *mapping_batch_new(array_list_t *fq_batch, pair_mng_t *pair_mng);
mapping_batch_t *mapping_batch_new_by_num(size_t num_reads, pair_mng_t *pair_mng);
void mapping_batch_free(mapping_batch_t *p);

//=====================================================================================
//=====================================================================================

/**
 * @brief Structure for store all process data in @a region_seeker_server.
 * 
 * Structure for store process data in @a region_seeker_server, contains an 
 * array list with all mappings found for all seeds in each read and a
 * batch with all reads unmapped.
 */
typedef struct cal_batch {
  array_list_t **allocate_mapping;  /**< array list with all mappings found for all seeds in each read*/
  fastq_batch_t *unmapped_batch;   /**< batch with all reads unmapped*/
} cal_batch_t;

//------------------------------------------------------------------------------------

/**
 * @brief  Constructor for the @a cal_batch_t structure.
 * @param  allocate_mapping_p array list with all mappings found for all seeds in each read
 * @param  unmapped_batch_p batch with all reads unmapped
 * @return Pointer to the new structure
 * 
 * Constructor function that allocates memory for
 * the cal_batch_t structure and initialize it.
 */
cal_batch_t* cal_batch_new(array_list_t **allocate_mapping, fastq_batch_t *unmapped_batch);

/**
 * @brief  Destrcutor for the @a cal_batch_t structure.
 * @param  allocate_cal_p[out] pointer to the structure to free
 * 
 * @a cal_batch_t destructor that frees the memory previosly
 * allocated by constructor @a cal_batch_new
 */
void cal_batch_free(cal_batch_t *allocate_cal);

//====================================================================================

/**
 * @brief Structure for store all process data in @a cal_seeker_server.
 * 
 * Structure for store process data in @a cal_seeker_server, store all cals found
 * for each read unmapped.
 */
typedef struct sw_batch {
  unsigned int num_reads;          /**< number of reads allocate in the batch*/
  array_list_t **allocate_cals_p;  /**< array list that store all cals found for each read */
  fastq_read_t **allocate_reads_p; /**< array that store for each read the header, the sequence and the quality*/
}sw_batch_t;

/**
 * @brief  Constructor for the @a sw_batch_t structure.
 * @param  num_reads number of reads that will be store
 * @param  allocate_cals_p array list that store all cals found for each read
 * @param  allocate_reads_p array that store for each read the header, the sequence and the quality
 * @return Pointer to the new structure
 * 
 * Constructor function that allocates memory for
 * the sw_batch_t structure and initialize it.
 */
sw_batch_t* sw_batch_new(unsigned int num_reads, array_list_t **allocate_cals_p, fastq_read_t **allocate_reads_p);

/**
 * @brief  Destrcutor for the @a sw_batch_t structure.
 * @param  sw_batch_p[out] pointer to the structure to free
 * 
 * @a sw_batch_t destructor that frees the memory previosly
 * allocated by constructor @a sw_batch_new
 */
void sw_batch_free(sw_batch_t *sw_batch_p);

//=====================================================================================

/**
 * @brief  Store in @a buffer_p all data information for one splice junction.
 * @param  chromosome chromosome where the splice junction was found
 * @param  strand strand where the splice junction was found
 * @param  start splice junction start coordinate 
 * @param  end splice junction end coordinate
 * @param  junction_id splice junction identification
 * @param  num_reads number of reads that support the splice junction
 * @param  buffer_p buffer to store splice junction data information  
 * @return Number of bytes stored in @a buffer_p
 * 
 * Store all data information for each splice junction found. It allocates all data 
 * in one buffer.
 */
unsigned int pack_junction(unsigned int chromosome, unsigned int strand, 
			   size_t start, size_t end, 
			   size_t junction_id, size_t num_reads, 
			   char *type, char* buffer_p);



typedef struct bwt_server_input bwt_server_input_t;
typedef struct region_seeker_input region_seeker_input_t;
typedef struct cal_seeker_input cal_seeker_input_t;
typedef struct preprocess_rna_input preprocess_rna_input_t;
typedef struct pair_server_input pair_server_input_t;
typedef struct sw_server_input sw_server_input_t;
typedef struct batch_writer_input batch_writer_input_t;


typedef struct batch {
  int mapping_mode;
  bwt_server_input_t *bwt_input;
  region_seeker_input_t *region_input;
  cal_seeker_input_t *cal_input;
  preprocess_rna_input_t *preprocess_rna;
  pair_server_input_t *pair_input;
  sw_server_input_t *sw_input;
  batch_writer_input_t *writer_input;
  mapping_batch_t *mapping_batch;
  void *optional_data;
} batch_t;


batch_t *batch_new(bwt_server_input_t *bwt_input,
                   region_seeker_input_t *region_input,
                   cal_seeker_input_t *cal_input,
                   pair_server_input_t *pair_input,
		   preprocess_rna_input_t *preprocess_rna,
                   sw_server_input_t *sw_input,
                   batch_writer_input_t *writer_input,
		   int mapping_mode,
                   mapping_batch_t *mapping_batch);

void batch_free(batch_t *b);

//======================================================================================

typedef struct buffer_item {
  fastq_read_t *read;
  array_list_t *items_list;
  void *aux_data;
} buffer_item_t;

buffer_item_t *buffer_item_new();
buffer_item_t *buffer_item_complete_new(fastq_read_t *read, array_list_t *cals_list, void *aux_data);
void buffer_item_insert_new_item(fastq_read_t *fq_read,
                                 linked_list_t *items_list,
                                 void *data,
                                 int type_items,
                                 linked_list_t *buffer,
                                 linked_list_t *buffer_hc,
                                 int phase);
void buffer_item_free(buffer_item_t *buffer_item);

void insert_file_item(fastq_read_t *fq_read, array_list_t *items, FILE *f_sa);

void insert_file_item_2(fastq_read_t *fq_read, array_list_t *items, FILE *f_hc);

//======================================================================================

//------------------------------------------------------------------

typedef struct sa_alignment {
  int left_close;
  int right_close;
  int num_sp;
  int complete;
  array_list_t *cals_list;
  cigar_code_t *c_left;
  cigar_code_t *c_right;
  cigar_code_t *c_final;
  int sp_middle[20];
  int reported;
} sa_alignment_t;

inline sa_alignment_t *sa_alignment_new(array_list_t *cals_list) {

  sa_alignment_t *sa_a = (sa_alignment_t *) malloc(sizeof(sa_alignment_t));

  sa_a->cals_list   = cals_list;
  sa_a->left_close  = 0;
  sa_a->right_close = 0;
  sa_a->c_left  = NULL;
  sa_a->c_right = NULL;

  memset(sa_a->sp_middle, 0, 20);
  sa_a->num_sp = 0;
  sa_a->complete = 0;
  sa_a->reported = 0;
  return sa_a;

}

inline void sa_alignment_free(sa_alignment_t *sa_alignment) {
  free(sa_alignment);
}

//------------------------------------------------------------------

typedef struct meta_alignment {
  int status;
  int sp_sw;
  int num_cigars;
  int type;
  int type_cigars[10];
  int score;
  float f_score;
  int flag;
  array_list_t *cals_list;
  cigar_code_t *middle_cigars[10];
  cigar_code_t *cigar_left;
  cigar_code_t *cigar_right;
  cigar_code_t *cigar_code;
} meta_alignment_t;

meta_alignment_t *meta_alignment_new();

typedef struct simple_alignment {
  int gap_start;
  int gap_end;
  int map_strand;
  int map_chromosome;
  size_t map_start;
  int map_distance;
  int cigar_len;
} simple_alignment_t;

typedef struct alignment_aux {
  int mapping_len;
  short int optional_fields_length;
  short int chromosome;
  int position;
  int map_quality;
  int num_cigar_operations;
  unsigned char seq_strand;	
  int cigar_len;
} alignment_aux_t;


void alignment_aux_init(alignment_t* alignment, alignment_aux_t *alignment_aux);


cal_t *convert_bwt_anchor_to_CAL(bwt_anchor_t *bwt_anchor, 
				 size_t read_start, size_t read_end);

void file_write_items(fastq_read_t *fq_read, array_list_t *items, 
		      unsigned char data_type, FILE *fd1, FILE *fd2, int mode);

fastq_read_t *file_read_fastq_reads(size_t *num_items, FILE *fd);

int file_read_cals(size_t num_items, array_list_t *list, 
		   fastq_read_t *fq_read, FILE *fd);

int file_read_meta_alignments(size_t num_items, array_list_t *list, 
			      fastq_read_t *fq_read, FILE *fd);

int file_read_alignments(size_t num_items, array_list_t *list, 
			 fastq_read_t *fq_read, FILE *fd);

int sa_file_read_alignments(size_t num_items, array_list_t *list, 
			    fastq_read_t *fq_read, FILE *fd);

void sa_file_write_alignments(fastq_read_t *fq_read, array_list_t *items, FILE *fd);

void sa_file_write_items(fastq_read_t *fq_read, array_list_t *items, FILE *fd);

#endif // BUFFERS_H
