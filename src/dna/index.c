#include "index.h"

//--------------------------------------------------------------------

index_t *index_new(char *dirname) {

    index_t *p = (index_t *) calloc(1, sizeof(index_t));

	char index_path[strlen(dirname) + 100];
	sprintf(index_path, "%s/params.txt", dirname);

	if (exists(index_path)) {
		// load SA index
		printf("Loading SA tables...\n");
		p->sa_index = sa_index3_new(dirname);
		p->mode = SA_MODE;
	} else {
		// load BWT index
		printf("Loading BWT tables...\n");
		p->bwt_index = bwt_index_new(dirname, false);
		p->genome = genome_new("dna_compression.bin", dirname, BWT_MODE);
		p->mode = BWT_MODE;
	}

	return p;
}

//--------------------------------------------------------------------

void index_free(index_t *p) {
    if (p) {
	    if (p->sa_index) sa_index3_free(p->sa_index);
	    if (p->bwt_index) bwt_index_free(p->bwt_index);
	    if (p->genome) genome_free(p->genome);
	    free(p);
	}
}

//--------------------------------------------------------------------

bam_header_t *index_create_bam_header(options_t *options, index_t *index) {
  bam_header_t *bam_header = NULL;
  if (index->mode == SA_INDEX) {
      bam_header = create_bam_header(options, index->sa_index->genome);
  } else {
      bam_header = create_bam_header_by_genome(index->genome);
  }
  return bam_header;
}

//--------------------------------------------------------------------

void write_sam_header_BWT(options_t *options, genome_t *genome, FILE *f);

void index_write_sam_header(options_t *options, FILE *file, index_t *index) {
  if (index->mode == SA_INDEX) {
      write_sam_header(options, index->sa_index->genome, file);
  } else {
      write_sam_header_BWT(options, index->genome, file);
  }
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

