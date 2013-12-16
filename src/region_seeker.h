#ifndef REGION_SEEKER_H
#define REGION_SEEKER_H

#include <stdio.h>

#include "containers/list.h"

#ifdef HPG_GPU
   #include "aligners/bwt/gpu.h"
#endif

#include "buffers.h"

#include "breakpoint.h"
#include "bs/methylation.h"

//====================================================================================
//  structures and prototypes
//====================================================================================

/**
 * @brief Structure for store all parameters needed for run @a region_seeker_server.
 * 
 * Structure for store some configuration values and data structures like lists 
 * to insert and read data batches.
 */
struct region_seeker_input{
  unsigned int region_threads;    /**< number of threads for run region seeker */ 
  unsigned int gpu_enable;
  list_t *unmapped_read_list_p;   /**< list for store batches with the reads no mapped */
  list_t* region_list_p;          /**< list to store batches with all regions found for each read */
  cal_optarg_t *cal_optarg_p;     /**< cal seeker configuration values */
  bwt_optarg_t *bwt_optarg_p;     /**< burrows wheeler transform configuration values */
  bwt_index_t *bwt_index_p;       /**< structure where were stored burrows wheeler transform index */
  bwt_index_t *bwt_index2_p;      /**< structure where were stored burrows wheeler transform second index */
  int padding_left;
  int padding_right;
  genome_t *genome;
  metaexons_t *metaexons;
  #ifdef HPG_GPU
    gpu_context_t *gpu_context;
  #endif
};

/**
 * @brief  Initializer for the @a region_seeker_input_t structure.
 * @param  unmapped_read_list_p list for store batches with the reads no mapped
 * @param  cal_optarg_p cal seeker configuration values 
 * @param  bwt_optarg_p burrows wheeler transform configuration values
 * @param  bwt_index_p structure where were stored burrows wheeler transform index
 * @param  region_list_p list to store batches with all regions found for each read
 * @param  region_threads number of threads for run region seeker
 *
 * 
 * Initialize all @a region_seeker_input_t fields with the input parameters.
 */
void region_seeker_input_init(list_t *unmapped_read_list_p, cal_optarg_t *cal_optarg_p, 
			      bwt_optarg_t *bwt_optarg_p, bwt_index_t *bwt_index_p, 
			      list_t* region_list_p, unsigned int region_threads, 
			      unsigned int gpu_enable, int padding_left, int padding_right, 
			      genome_t *genome, metaexons_t *metaexons,
			      region_seeker_input_t *input_p);
//--------------------------------------------------------------------------------------

/**
 * @brief Region Seeker server generate regions. 
 * @param input_p all configuration values and data structures needed by the server
 *
 * Region Seeker server extract batches of reads from @a unmapped_read_list_p and
 * for each read call @a bwt_map_exact_seeds_seq for make seeds and search mappings
 * of these. Finally all reads with their mappings are stored in @a region_list_p 
 * for continue process these in next phases of pipeline.      
 */
int region_seeker_server(region_seeker_input_t *input);
//====================================================================================

int apply_seeding(region_seeker_input_t* input, batch_t *batch);
int apply_seeding_bs(region_seeker_input_t* input, batch_t *batch);

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

#endif // REGION_SEEKER_H

