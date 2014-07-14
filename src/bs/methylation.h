#ifndef METHYLATION_H
#define METHYLATION_H

//====================================================================================

#include <stdio.h>
#include <stdlib.h>

#include "commons/commons.h"
#include "commons/system_utils.h"

#include "containers/array_list.h"
#include "containers/list.h"

#include "aligners/bwt/bwt.h"
#include "aligners/bwt/genome.h"
#include "aligners/sw/smith_waterman.h"

#include "bioformats/fastq/fastq_file.h"

#include "buffers.h"
#include "pair_server.h"

#include "hash_table.h"
#include "array_list_bs.h"

//====================================================================================

// values for the diferent combination of nucleotides
#define ACGT 0
#define AGT  1
#define ACT  2
#define AT   3

// limites to use with the histogram data of the reads
#define LIMIT_INF   0.20
#define LIMIT_SUP   0.30


// macros for postproces

#define search_gen()							\
  ({res = tmp & (unsigned long long)3;					\
    switch(res){							\
    case 0:								\
      break;								\
    case 1:								\
      if (read[pos] == c1){						\
	postprocess_bs(query_name, '+', c, init + pos, 'Z', strand, 0, bs_context->context_CpG); \
	bs_context->CpG_methyl++;					\
      } else {								\
	if (read[pos] == c2) {					\
	  postprocess_bs(query_name, '-', c, init + pos, 'z', strand, 0, bs_context->context_CpG); \
	  bs_context->CpG_unmethyl++;					\
	} else {							\
	  postprocess_bs(query_name, '.', c, init + pos, 'M', strand, 3, bs_context->context_MUT); \
	  bs_context->MUT_methyl++;					\
	}								\
      }									\
      break;								\
    case 2:								\
      if (read[pos] == c1) {					\
	postprocess_bs(query_name, '+', c, init + pos, 'X', strand, 0, bs_context->context_CHG); \
	bs_context->CHG_methyl++;					\
      } else {								\
	if (read[pos] == c2) {					\
	  postprocess_bs(query_name, '-', c, init + pos, 'x', strand, 0, bs_context->context_CHG); \
	  bs_context->CHG_unmethyl++;					\
	} else {							\
	  postprocess_bs(query_name, '.', c, init + pos, 'M', strand, 3, bs_context->context_MUT); \
	  bs_context->MUT_methyl++;					\
	}								\
      }									\
      break;								\
    case 3:								\
      if (read[pos] == c1) {						\
	postprocess_bs(query_name, '+', c, init + pos, 'H', strand, 2, bs_context->context_CHH); \
	bs_context->CHH_methyl++;					\
      } else {								\
	if (read[pos] == c2) {					\
	  postprocess_bs(query_name, '-', c, init + pos, 'h', strand, 2, bs_context->context_CHH); \
	  bs_context->CHH_unmethyl++;					\
	} else {							\
	  postprocess_bs(query_name, '.', c, init + pos, 'M', strand, 3, bs_context->context_MUT); \
	  bs_context->MUT_methyl++;					\
	}								\
      }									\
      break;								\
    }									\
    tmp = tmp >> 2;})


//====================================================================================

typedef struct metil_file {
  char* filenameCpG;                     /**< Metilation file name for CpG regions. */
  char* filenameCHG;                     /**< Metilation file name for CpG regions. */
  char* filenameCHH;                     /**< Metilation file name for CpG regions. */
  char* filenameMUT;                     /**< Metilation file name for Mutations.   */
  char* filenameSTAT;                    /**< Metilation file name for statistics.  */
  char* mode;                            /**< Open mode ("r", "w").                 */ // not in use
  genome_t *genome;                      /**< Reference to the original genome.     */

  FILE *CpG;                             /**< File pointer to CpG output.           */
  FILE *CHG;                             /**< File pointer to CHG output.           */
  FILE *CHH;                             /**< File pointer to CHH output.           */
  FILE *MUT;                             /**< File pointer to Mutation output.      */
  FILE *STAT;                            /**< File pointer to Statistics output.    */

  hash_table_t *table_isles;             /**< Structure with the hash table values  */

  // check the values in order to change the type of the counter variables
  size_t CpG_methyl;                 /**< Global Counter for methylated Cytosines in CpG context   */
  size_t CpG_unmethyl;               /**< Global Counter for unmethylated Cytosines in CpG context */
  size_t CHG_methyl;                 /**< Global Counter for methylated Cytosines in CHG context   */
  size_t CHG_unmethyl;               /**< Global Counter for unmethylated Cytosines in CHG context */
  size_t CHH_methyl;                 /**< Global Counter for methylated Cytosines in CHH context   */
  size_t CHH_unmethyl;               /**< Global Counter for unmethylated Cytosines in CHH context */
  size_t MUT_methyl;                 /**< Global Counter for mutated Cytosines                     */
  size_t num_bases;                  /**< Global Counter for number of bases in the batch          */
} metil_file_t;

//====================================================================================

/**
 * @brief  Initialization of the metilation structure.
 * @param  metil_file metilation structure to use in the final analysis.
 * @param  dir        output directory.
 * @param  genome     pointer to the compressed original genome.
 * 
 * Initialize the variables in the structure, and open the output files.
 */
void metil_file_init(metil_file_t *metil_file, char *dir, genome_t *genome);

/**
 * @brief  Release the memory.
 * @param  metil_file metilation structure.
 * 
 * Release the memory used by the structure, and close the files in use.
 */
void metil_file_free(metil_file_t *metil_file);

//====================================================================================

// struct ported to the array_list_bs file to avoid conflicts

/*
typedef struct metil_data {
  char*  query_name;                    //
  char   status;                        //
  int    chromosome;                    //
  size_t start;                         //
  char   context;                       //
  int    strand;                        //
  int    zone;                          //
} metil_data_t;
*/

//====================================================================================

/**
 * @brief  
 * @param  
 * @param  
 * @param  
 * 
 * 
 */
void metil_data_init(metil_data_t *metil_data, char *query, char status, int chromosome, size_t start, char context, int strand, int zone);

/**
 * @brief  
 * @param  
 * 
 * 
 */
void metil_data_free(metil_data_t *metil_data);

//====================================================================================

/**
 * @brief  Transfor the sequence to the metilation case
 * @param  refs sequence of reference
 * @param  len length of the sequence
 * @param  type type of transformation
 * 
 * Replace the bases of the sequence depending the type.
 */
void replace(char *refs, int len, int type);

//====================================================================================

/**
 * @brief  Obtain the complementary.
 * @param  c base of a sequence.
 * @return complementary base.
 * 
 * Obtain the complementary base of the input.
 */
char complement (char c);

//====================================================================================

/**
 * @brief  Obtain the reverse complementary of a dna sequence.
 * @param  orig sequence to transform.
 * @param  dest sequence transformed.
 * @param  len  length of the sequence.
 * 
 * Transform the @a orig sequence of DNA, to the reverse complementary.
 */
void rev_comp(char *query, char *seq, int len);

//====================================================================================

/**
 * @brief  Transform the sequence into the complementary.
 * @param  seq sequence of nucleotides.
 * @param  len length of the sequence.
 * 
 * Complement the sequence element by element.
 */
void comp(char *seq, int len);

//====================================================================================

/**
 * @brief  Make a copy of an array_list.
 * @param  dest pointer to the new structure that is gona store the reads.
 * @param  src  pointer to the structure that is gone be copied.
 * 
 * Copy an array of reads into a new one.
 */
void copy_array(array_list_t *dest, array_list_t *src);

//====================================================================================

/**
 * @brief  Make four copies of the array of reads for the bisulfite case.
 * @param  src   original data array.
 * @param  dest1 destination array 1.
 * @param  dest2 destination array 2.
 * @param  dest3 destination array 3.
 * @param  dest4 destination array 4.
 * 
 * Copy the reads and insert in the four arrays (old version, no use)
 */
void cpy_array_bs(array_list_t *src, array_list_t *dest1, array_list_t *dest2, array_list_t *dest3, array_list_t *dest4);

//====================================================================================

/**
 * @brief  Make four copies of the array of reads, and transform for the bisulfite case.
 * @param  src original data array.
 * @param  dest_ct     destination to the ct transformation.
 * @param  dest_ct_rev destination to the ct rev comp transformation.
 * @param  dest_ga     destination to the ga transformation.
 * @param  dest_ga_rev destination to the ga rev comp transformation.
 * 
 * Copy the reads and transform it before insert in the new array.
 */
void cpy_transform_array_bs(array_list_t *src, array_list_t *dest_ct, array_list_t *dest_ct_rev,
			    array_list_t *dest_ga, array_list_t *dest_ga_rev);

//====================================================================================

/**
 * @brief  Combine the alignments from both arrays.
 * @param  dest array is gona store the combination.
 * @param  src array is going to unify the other.
 * 
 * Insert the alignments from @a src into the mappings on @a dest.
 */
void insert_mappings_array(array_list_t **dest, array_list_t **src);

//====================================================================================

/**
 * @brief  Include the new mappings existing in the "src" into "dest".
 * @param  dest Mappings that are gone store the combination.
 * @param  src mapping are goin to be joined.
 * 
 * Unify two mapping_list into one.
 */
void insert_mappings(array_list_t *dest, array_list_t *src);

//====================================================================================

/**
 * @brief  Include the new regions existing in the "src" into "dest".
 * @param  dest Regions that are gone store the combination.
 * @param  src regions are goin to be joined.
 * 
 * Unify two region_list into one.
 */
void insert_regions(array_list_t *dest, array_list_t *src);

//====================================================================================

/**
 * @brief  Set the sequence strand of the mappings on the array to 1.
 * @param  src mappings to transform.
 * 
 * Set all the mappins from all the reads in @a src on the reverse strand.
 */
void transform_mappings_array(array_list_t **src);

//====================================================================================

/**
 * @brief  Set the sequence strand of the mappings to 1.
 * @param  src Array of mappings.
 * 
 * Set all the mappins in @a src on the reverse strand.
 */
void transform_mappings(array_list_t *src);

//====================================================================================

/**
 * @brief  Set the region strand to 1.
 * @param  src Array of regions.
 * 
 * Set all the regions in @a src on the reverse strand.
 */
void transform_regions(array_list_t *src);

//====================================================================================

/**
 * // not used
 * @brief  Select the index values where no read is mapped.
 * @param  num_unmapped  store the number of unmapped reads.
 * @param  unmapped_indices store the index of unmapped reads.
 * @param  indices_1 vector with the count of unmapped reads.
 * @param  num_reads number of reads analyzed.
 * @param  lists,lists2 arrays with the information relative to the mappings.
 * 
 * Updates the unmapped index for the next stage.
 */
void select_targets_bs(size_t *num_unmapped, size_t *unmapped_indices,
		       size_t *indices_1, size_t num_reads,
		       array_list_t **lists, array_list_t **lists2);

//====================================================================================

/**
 * // not used
 * @brief  
 * @param  num_reads        
 * @param  unmapped_indices 
 * @param  unmapped         
 * @param  indices1         
 * 
 * 
 */
void update_targets_bs(size_t num_reads, size_t *unmapped_indices,
		       size_t unmapped, size_t *indices);

//====================================================================================

/**
 * @brief  Set the original read in the sequences of each mapping.
 * @param  src1  Array with the mappings of the reads transformed.
 * @param  src1  Array with the mappings of the reads transformed.
 * @param  orig Fastq batch with the original reads.
 * 
 * Go on the array of mappings, and set the sequence to the original.
 */
void revert_mappings_seqs(array_list_t **src1, array_list_t **src2, array_list_t *orig);

//====================================================================================

/**
 * @brief  Obtain the reverse complementary reads stored in the array
 * @param  dest 
 * @param  src 
 * 
 * 
 */
void rev_comp_array(array_list_t *dest, array_list_t *src);

//====================================================================================

/**
 * @brief  Make an histogram of the read
 * @param  seq  
 * @param  len  
 * @param  type 
 * 
 * Return true if the number of certain base is greater than the filter
 */
int histogram_seq(char *seq, size_t len, int type);

//====================================================================================

/**
 * @brief  Write the metilation status of each Cytosine
 * @param  array_list 
 * @param  metil_file 
 * 
 * 
 */
void write_metilation_status(array_list_t *array_list, metil_file_t *metil_file);

//====================================================================================

/*
 * @brief  Write the metilation status of each Cytosine
 * @param  array_list array with the information to report
 * @param  metil_file metilation structura with the output information
 * 
 * Report in diferent files depending of the region, the status of the Cytosines aligned in the pipeline
 */
//void methylation_status_report(metil_file_t *metil_file, batch_t *batch);

//====================================================================================

/**
 * @brief  Stage to determine the status of each Cytosine in the sequence alignment                                       
 * @param  batch                                                                                                           
 *                                                                                                                         
 *                                                                                                                         
 */
int methylation_status_report(sw_server_input_t* input, batch_t *batch);

//====================================================================================                                     

/**
 * @brief  Add the metilation status of each Cytosine to the write list                                                    
 * @param  array_list                                                                                                      
 * @param  bs_status                                                                                                       
 *                                                                                                                         
 *                                                                                                                         
 */
//void add_metilation_status(array_list_t *array_list, array_list_t *list, bs_context_t *bs_context, genome_t * genome, array_list_t * orig_seq, size_t index, int conversion);
void add_metilation_status(array_list_t *array_list, bs_context_t *bs_context, genome_t * genome, array_list_t * orig_seq, size_t index, int conversion);

//====================================================================================                                     

/**
 * @brief  Add the metilation status of each Cytosine to the write list
 * @param  array_list
 * @param  bs_status
 *
 *
 */
void add_metilation_status_bin(array_list_t *array_list, bs_context_t *bs_context,
			       unsigned long long **gen_binCT, unsigned long long **gen_binGA,
			       array_list_t * orig_seq, size_t index, int conversion);

//====================================================================================                                     

/**
 * @brief  
 * @param  
 * @param  
 *
 *
 */
void search_methylation(int c, size_t init, size_t end, unsigned long long **values, char *read, bs_context_t *bs_context, char strand, int type, char *query_name);

//====================================================================================                                     

/**
 * @brief  
 * @param  
 * @param  
 *                                                                                                                         
 *                                                                                                                         
 */
void postprocess_bs(char *query_name, char status, size_t chromosome, size_t start, char context, int strand, int region,
		    array_list_t *list);

//====================================================================================                                     

/**
 * @brief  
 * @param  
 * @param  
 *                                                                                                                         
 *                                                                                                                         
 */
void postproc_bs(char *query_name, char status, size_t chromosome, size_t start, char context, int strand, int region,
		 array_list_bs_t *list);

//====================================================================================                                     

/**
 * @brief  Write the metilation status of each Cytosine
 * @param  array_list 
 * @param  metil_file 
 * 
 * 
 */
void write_bs_context(metil_file_t *metil_file, bs_context_t *bs_context);

//====================================================================================

int encode_context(char* filename, char* directory);

//====================================================================================

int load_encode_context(char* directory, unsigned long long **valuesCT, unsigned long long **valuesGA);

//====================================================================================

#endif // METHYLATION_H
