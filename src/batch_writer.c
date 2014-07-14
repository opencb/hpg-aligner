#include "batch_writer.h"

//------------------------------------------------------------------------------------

bam_header_t *create_bam_header_by_genome(genome_t *genome) {

  bam_header_t *bam_header = (bam_header_t *) calloc(1, sizeof(bam_header_t));

  int num_targets = genome->num_chromosomes;

  bam_header->n_targets = num_targets;
  bam_header->target_name = (char **) calloc(num_targets, sizeof(char *));
  bam_header->target_len = (uint32_t*) calloc(num_targets, sizeof(uint32_t));
  bam_header->text = strdup("@PG\tID:HPG-Aligner\tVN:1.0\n");
  for (int i = 0; i < num_targets; i++) {
    bam_header->target_name[i] = strdup(genome->chr_name[i]);
    bam_header->target_len[i] = genome->chr_size[i] + 1;
  }
  bam_header->l_text = strlen(bam_header->text);

  return bam_header;
}

//------------------------------------------------------------------------------------

void batch_writer_input_init(char* match_filename, char* splice_exact_filename, 
			     char* splice_extend_filename, 
			     linked_list_t* list_p, genome_t* genome, 
			     batch_writer_input_t* input_p) {

  input_p->match_filename = match_filename;
  input_p->splice_exact_filename = splice_exact_filename;
  input_p->splice_extend_filename = splice_extend_filename;
  input_p->list_p = list_p;
  input_p->genome = genome;

  // internal
  input_p->bam_file = NULL;
  input_p->total_batches = 0;
  input_p->total_reads = 0;
  input_p->total_mappings = 0;
  input_p->num_mapped_reads = 0;
  input_p->limit_print = 10000;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
