#include "methylation.h"

const size_t margin_bwt = 10;
const size_t margin_sw  = 25;

//====================================================================================

void replace(char * refs, int len, int type) {
  char c1, c2;

  switch (type) {
  case ACGT:
    // case with no transformation required
    return;
  case AGT:
    c1 = 'C';
    c2 = 'T';
    break;
  case ACT:
    c1 = 'G';
    c2 = 'A';
    break;
  case AT:
    // apply both transformations, and end the execution
    replace(refs, len, 1);
    replace(refs, len, 2);
    return;
  default:
    // if the type value is not recognised, nothing is done
    return;
  }

  //printf("replace\nleng = %lu\ntype = %lu\n", len, type);
  //printf("c1 = %c\nc2 = %c\n\n", c1, c2);

  // transforms the refs sequence, replacing the caracter c1 with the caracter c2
  for (int j = 0; j < len; j++) {
    if (refs[j] == c1) {
      refs[j] = c2;
    }
  }

  return;
}

//====================================================================================

char complement (char c) {
  // return the complementary base
  switch (c) {
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  default:
    return c;
  }
}

//====================================================================================

void rev_comp(char *orig, char *dest, int len) {
  // put the reverse complementary sequence of orig in dest
  for (int i = 0; i < len; i++) {
    dest[len - i - 1] = complement(orig[i]);
  }
  dest[len] = '\0';
}

//====================================================================================

void comp(char *seq, int len) {
  // obtain the complementary sequence of seq
  for (int i = 0; i < len - 1; i++) {
    seq[i] = complement(seq[i]);
  }
}

//====================================================================================

void replace_array(array_list_t *reads, int type) {
  size_t num_reads = array_list_size(reads);
  fastq_read_t* fq_read;

  //printf("reads = %lu\n", num_reads);
  
  // replace the reads in the array depending the value of type
  for (size_t i = 0; i < num_reads; i++) {
    fq_read = (fastq_read_t *) array_list_get(i, reads);

    //printf("read = %lu\n", i);
    //printf("src = %s\n", fq_read->sequence);
    replace(fq_read->sequence, fq_read->length, type);
    //printf("src = %s\ntam = %lu\ntype = %lu\n\n", fq_read->sequence, fq_read->length, type);
  }
}

//====================================================================================

void rev_comp_array(array_list_t *dest, array_list_t *src) {
  size_t num_reads = array_list_size(src);
  fastq_read_t* fq_read_src;
  fastq_read_t* fq_read_dest;

  // get the reverse complementary sequence in all the array
  for (size_t i = 0; i < num_reads; i++) {
    fq_read_src  = (fastq_read_t *) array_list_get(i, src);
    fq_read_dest = (fastq_read_t *) array_list_get(i, dest);

    //printf("read = %lu\n", i);
    //printf("src = %s\ndst = %s\ntam = %lu\n\n", fq_read_src->sequence, fq_read_dest->sequence, fq_read_src->length);
    rev_comp(fq_read_src->sequence, fq_read_dest->sequence, fq_read_src->length);
    //printf("src = %s\ndst = %s\ntam = %lu\n\n", fq_read_src->sequence, fq_read_dest->sequence, fq_read_src->length);
  }
}

//====================================================================================

void copy_array(array_list_t *dest, array_list_t *src) {
  size_t num_reads = array_list_size(src);
  fastq_read_t* fq_read_src;
  fastq_read_t* fq_read_dest;
  
  for (size_t i = 0; i < num_reads; i++) {
    fq_read_src  = (fastq_read_t *) array_list_get(i, src);
    //fq_read_dest = fastq_read_new(fq_read_src);
    //insert_in_array(dest, fq_read_dest);
  }
}

//====================================================================================

void cpy_array_bs(array_list_t *src, array_list_t *dest1, array_list_t *dest2, array_list_t *dest3, array_list_t *dest4) {
  size_t num_reads = array_list_size(src);
  fastq_read_t* fq_read_src;
  fastq_read_t* fq_read_dest;
  
  // make four copies of the original array
  for (size_t i = 0; i < num_reads; i++) {
    fq_read_src  = (fastq_read_t *) array_list_get(i, src);
    fq_read_dest = fastq_read_dup(fq_read_src);
    array_list_insert(fq_read_dest, dest1);
    fq_read_dest = fastq_read_dup(fq_read_src);
    array_list_insert(fq_read_dest, dest2);
    fq_read_dest = fastq_read_dup(fq_read_src);
    array_list_insert(fq_read_dest, dest3);
    fq_read_dest = fastq_read_dup(fq_read_src);
    array_list_insert(fq_read_dest, dest4);
  }
}

//====================================================================================

void cpy_transform_array_bs(array_list_t *src, array_list_t *dest_ct, array_list_t *dest_ct_rev, array_list_t *dest_ga, array_list_t *dest_ga_rev) {
  size_t num_reads = array_list_size(src);
  fastq_read_t *fq_read_src;
  fastq_read_t *fq_read_dest;
  fastq_read_t *fq_read_tmp;

  // read element by element in the array, transform each one, and put in the new arrays
  for (size_t i = 0; i < num_reads; i++) {
    fq_read_src  = (fastq_read_t *) array_list_get(i, src);

    fq_read_dest = fastq_read_dup(fq_read_src);
    replace(fq_read_dest->sequence, fq_read_dest->length, AGT);
    array_list_insert(fq_read_dest, dest_ct);

    fq_read_tmp = fastq_read_dup(fq_read_src);
    rev_comp(fq_read_dest->sequence, fq_read_tmp->sequence, fq_read_dest->length);
    array_list_insert(fq_read_tmp, dest_ct_rev);

    fq_read_dest = fastq_read_dup(fq_read_src);
    replace(fq_read_dest->sequence, fq_read_dest->length, ACT);
    array_list_insert(fq_read_dest, dest_ga);

    fq_read_tmp = fastq_read_dup(fq_read_src);
    rev_comp(fq_read_dest->sequence, fq_read_tmp->sequence, fq_read_dest->length);
    array_list_insert(fq_read_tmp, dest_ga_rev);
  }
}

//====================================================================================

void insert_mappings_array(array_list_t **dest, array_list_t **src) {
  size_t num_reads = array_list_size(src);
  size_t num_mappings;
  alignment_t *align_tmp;

  for (size_t i = 0; i < num_reads; i++) {
    //isnert_mappings(dest[i], src[i]);
    num_mappings = array_list_size(src[i]);
    for (size_t j = 0; j < num_mappings; j++) {
      align_tmp  = (alignment_t *) array_list_get(j, src[i]);
      // insert the alignments from the 'src' into 'dest'
      array_list_insert(align_tmp, dest[i]);
    }
  }
}

//====================================================================================

void insert_mappings(array_list_t *dest, array_list_t *src) {
  size_t num_mappings;
  alignment_t *align_tmp;

  num_mappings = array_list_size(src);
  for (size_t j = 0; j < num_mappings; j++) {
    align_tmp  = (alignment_t *) array_list_get(j, src);
    // insert the alignments from the 'src' into the 'dest'
    array_list_insert(align_tmp, dest);
  }
}

//====================================================================================

void insert_regions(array_list_t *dest, array_list_t *src) {
  size_t num_mappings;
  region_t *region_tmp;

  num_mappings = array_list_size(src);
  for (size_t j = 0; j < num_mappings; j++) {
    region_tmp  = (region_t *) array_list_get(j, src);
    // insert the alignments from the 'src' into the 'dest'
    array_list_insert(region_tmp, dest);
  }
}

//====================================================================================

void transform_mappings_array(array_list_t **src){
  size_t num_reads = array_list_size(src);
  size_t num_mappings;
  alignment_t *align_tmp;

  // go over all the sequences
  for (size_t i = 0; i < num_reads; i++) {
    num_mappings = array_list_size(src[i]);
    // go over all the alignments
    for (size_t j = 0; j < num_mappings; j++) {
      align_tmp  = (alignment_t *) array_list_get(j, src[i]);
      // modify the strand, and make it reverse
      align_tmp->seq_strand = 1;
    }
  }
}

//====================================================================================

void transform_mappings(array_list_t *src){
  size_t num_mappings;
  alignment_t *align_tmp;

  num_mappings = array_list_size(src);
  // go over all the alignments
  for (size_t j = 0; j < num_mappings; j++) {
    align_tmp  = (alignment_t *) array_list_get(j, src);
    // modify the strand, and make it reverse
    align_tmp->seq_strand = 1;
  }
}

//====================================================================================
void transform_regions(array_list_t *src){
  size_t num_mappings;
  region_t *region_tmp;

  num_mappings = array_list_size(src);
  // go over all the regions
  for (size_t j = 0; j < num_mappings; j++) {
    region_tmp  = (region_t *) array_list_get(j, src);
    // modify the strand, and make it reverse
    region_tmp->strand = 1;
  }
}

//====================================================================================

void select_targets_bs(size_t *num_unmapped, size_t *unmapped_indices,
		       size_t *indices_1, size_t num_reads,
		       array_list_t **lists, array_list_t **lists2) {

  *num_unmapped = 0;

  for (size_t i = 0; i < num_reads; i++) {
    if (indices_1[i] == 4) {
      unmapped_indices[(*num_unmapped)++] = indices_1[i];
      array_list_set_flag(0, lists[indices_1[i]]);
      array_list_set_flag(0, lists2[indices_1[i]]);
    }
  }
}

//====================================================================================

void update_targets_bs(size_t num_reads, size_t *unmapped_indices,
		       size_t unmapped, size_t *indices) {
  for (size_t i = 0; i < unmapped; i++) {
    if (indices[i] < num_reads) {
      unmapped_indices[indices[i]]++;
    }
  }
}

//====================================================================================

void revert_mappings_seqs(array_list_t **src1, array_list_t **src2, array_list_t *orig) {
  //printf("+++++++++++++++ is orig NULL ? %i, size %i\n", (orig == NULL), array_list_size(orig));

  size_t num_mappings;
  alignment_t *align_tmp;
  fastq_read_t *fastq_orig;
  size_t num_reads = array_list_size(orig);

  //printf("num reads = %lu\t", num_reads);

  //printf("num reads = %lu\t", array_list_size(orig));
  //printf("num reads src1 = %lu\t", array_list_size(src1));
  //printf("num reads src2 = %lu\t", array_list_size(src2));

  // go over all the sequences
  for (size_t i = 0; i < num_reads; i++) {
    fastq_orig  = (fastq_read_t *) array_list_get(i, orig);
 
    // go over all the alignments in list 1
    num_mappings = array_list_size(src1[i]);
    for (size_t j = 0; j < num_mappings; j++) {
      align_tmp  = (alignment_t *) array_list_get(j, src1[i]);
      // free existing memory and copy the original
      if (align_tmp->sequence != NULL) free(align_tmp->sequence);

      align_tmp->sequence = strdup(fastq_orig->sequence);
    }

    // go over all the alignments in list 2
    num_mappings = array_list_size(src2[i]);
    for (size_t j = 0; j < num_mappings; j++) {
      align_tmp  = (alignment_t *) array_list_get(j, src2[i]);
      // free existing memory and copy the original
      if (align_tmp->sequence != NULL) free(align_tmp->sequence);

      align_tmp->sequence = strdup(fastq_orig->sequence);
    }
  }
}

//====================================================================================

char *obtain_seq(alignment_t *alig, fastq_read_t * orig) {
  //char *read = alig->sequence;
  char *read = orig->sequence;
  char *cigar = strdup(alig->cigar);
  //char *cigar = strdup("100M1D");
  int num;
  char car;
  int cont, pos, pos_read;
  int operations;
  int len = strlen(cigar) - 1;
  char *seq = (char *)calloc(1024, sizeof(char));
  //printf("cigar %s\n", cigar);
  //printf("read  %s\n", read);
  
  pos = 0;
  pos_read = 0;
  
  for (operations = 0; operations < alig->num_cigar_operations; operations++) {
    sscanf(cigar, "%i%c%s", &num, &car, cigar);
    //printf("%3i %c %s\n",  num, car, cigar);
    if (car == 'M' || car == '=' || car == 'X') {
      for (cont = 0; cont < num; cont++, pos++, pos_read++) {
	seq[pos] = read[pos_read];
      }
    }
    else {
      if (car == 'D' || car == 'N') {
	pos_read += num - 1;
      }
      else {
	if (car == 'I' || car == 'H' || car == 'S') {
	  for (cont = 0; cont < num; cont++, pos++) {
	    seq[pos] = '-';
	  }
	}
      }
    }
  }
  seq[pos] = '\0';

  //printf("seq   %s\n", seq);

  free(cigar);
  return seq;
}

//====================================================================================

void write_metilation_status_new(array_list_t *array_list, metil_file_t *metil_file) {

  //printf("Init metilation status\n");
  /*
  printf("CpG %s\nCHG %s\nCHH %s\nMUT %s\n\n",
	 metil_file->filenameCpG, metil_file->filenameCHG, metil_file->filenameCHH, metil_file->filenameMUT);
  */  

  FILE * CpG = metil_file->CpG;
  if (CpG == NULL) {
    printf("reopen CpG file\n");
    CpG =fopen(metil_file->filenameCpG, "a");
  }
  FILE * CHG = metil_file->CHG;
  if (CHG == NULL) {
    printf("reopen CHG file\n");
    CHG = fopen(metil_file->filenameCHG, "a");
  }
  FILE * CHH = metil_file->CHH;
  if (CHH == NULL) {
    printf("reopen CHH file\n");
    CHH = fopen(metil_file->filenameCHH, "a");
  }
  FILE * MUT = metil_file->MUT;
  if (MUT == NULL) {
    printf("reopen CHH file\n");
    MUT = fopen(metil_file->filenameMUT, "a");
  }

  size_t num_items = array_list_size(array_list);
  alignment_t *alig;
  genome_t *genome = metil_file->genome;
  hash_table_t *table_bs = metil_file->table_isles;
  char *seq, *gen, *read;
  size_t len, end, start;

  char *cigar;
  int num;
  char car;
  int cont, pos, pos_read;

  int file_error;

  for (size_t j = 0; j < num_items; j++) {
    alig = (alignment_t *) array_list_get(j, array_list);
    if (alig != NULL && alig->is_seq_mapped) {
      /*
      printf("alignment %lu\n", j);
      printf("query_name %s\n", alig->query_name);
      printf("sequence   %s\n", alig->sequence);
      printf("cigar      %s\n", alig->cigar);
      printf("position   %i\n", alig->position);
      printf("mate pos   %i\n", alig->mate_position);
      printf("temp len   %i\n", alig->template_length);
      printf("chromo     %i\n", alig->chromosome);
      printf("strand     %i\n", alig->seq_strand);
      printf("mate  chro %i\n", alig->mate_chromosome);
      printf("map qual   %i\n", alig->map_quality);
      printf("n cigar    %i\n", alig->num_cigar_operations);
      */      

      /*      
      seq = obtain_seq(alig);
      len = strlen(seq);
      */
      len = strlen(alig->sequence);
      gen = (char *)calloc(len + len + 2, sizeof(char));

      start = alig->position + 1;
      end = start + len - 1;
      if (end >= genome->chr_size[alig->chromosome]) {
        end = genome->chr_size[alig->chromosome] - 1;
      }

      genome_read_sequence_by_chr_index(gen, alig->seq_strand, alig->chromosome, &start, &end, genome);

      //printf("\nseq %s\ngen %s\n", seq, gen);

      read = alig->sequence;
      cigar = strdup(alig->cigar);
      for (int operations = 0; operations < alig->num_cigar_operations; operations++) {
	sscanf(cigar, "%i%c%s", &num, &car, cigar);
	//printf("%3i %c %s\n",  num, car, cigar);
	if (car == 'M' || car == '=') {
	  for (cont = 0; cont < num; cont++, pos++, pos_read++) {
	    //seq[pos] = read[pos_read];
	    if (read[pos_read] == 'C') {
	      if (gen[pos] == 'C') {
		// methylated cytosine
		if (gen[pos + 1] == 'G') {
		  // CpG zone
		  //printf("%s\t+\t%i %i\t%lu\tZ\n", alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
		  /*
		  file_error = fprintf(CpG, "%s\t+\t%i %i\t%lu\tZ\n",
				       alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
		  if (file_error < 0) {
		    printf("Error al escribir\n");
		    exit(-1);
		  }
		  */
		}
		else {
		  if (gen[pos + 2] == 'G') {
		    // CHG zone
		    //printf("%s\t+\t%i %i\t%lu\tX\n", alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
		    /*
		    file_error = fprintf(CHG, "%s\t+\t%i %i\t%lu\tX\n",
					 alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
		    if (file_error < 0) {
		      printf("Error al escribir\n");
		      exit(-1);
		    }
		    */
		  }
		  else {
		    // CHH zone
		    //printf("%s\t+\t%i %i\t%lu\tH\n", alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
		    /*
		    file_error = fprintf(CHH, "%s\t+\t%i %i\t%lu\tH\n",
					 alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
		    if (file_error < 0) {
		      printf("Error al escribir\n");
		      exit(-1);
		    }
		    */
		  }
		}
	      }
	      else {
		// mutated cytosine
		//printf("%s\t+\t%i %i\t%lu\tM\n", alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
		/*
		file_error = fprintf(MUT, "%s\t+\t%i %i\t%lu\tM\n",
				     alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
		if (file_error < 0) {
		  printf("Error al escribir\n");
		  exit(-1);
		}
		*/
	      }
	    }
	  }
	}
	else {
	  if (car == 'D') {
	    pos_read += num - 1;
	  }
	  else {
	    if (car == 'I') {
	      for (cont = 0; cont < num; cont++, pos++) {
		//seq[pos] = 'N';
	      }
	    }
	  }
	}
      }
      free(cigar);

      if (seq) free(seq);
      if (gen) free(gen);
    }
  }
  //printf("end status\n");
}

//====================================================================================

//void metil_file_init(metil_file_t *metil_file, char *CpG, char *CHG, char *CHH, char *MUT, genome_t *genome) {
void metil_file_init(metil_file_t *metil_file, char *dir, genome_t *genome) {
  char *name_tmp = malloc(128 * sizeof(char));

  int file_error;

  sprintf(name_tmp, "%s/CpG_context.txt", dir);
  metil_file->filenameCpG = strdup(name_tmp);
  sprintf(name_tmp, "%s/CHG_context.txt", dir);
  metil_file->filenameCHG = strdup(name_tmp);
  sprintf(name_tmp, "%s/CHH_context.txt", dir);
  metil_file->filenameCHH = strdup(name_tmp);
  sprintf(name_tmp, "%s/MUT_context.txt", dir);
  metil_file->filenameMUT = strdup(name_tmp);
  sprintf(name_tmp, "%s/Statistics.txt", dir);
  metil_file->filenameSTAT = strdup(name_tmp);

  /*
  printf("CpG %s\nCHG %s\nCHH %s\nMUT %s\n\n",
	 metil_file->filenameCpG, metil_file->filenameCHG, metil_file->filenameCHH, metil_file->filenameMUT);
  */

  metil_file->genome = genome;

  metil_file->CpG  = fopen(metil_file->filenameCpG,  "w");
  metil_file->CHG  = fopen(metil_file->filenameCHG,  "w");
  metil_file->CHH  = fopen(metil_file->filenameCHH,  "w");
  metil_file->MUT  = fopen(metil_file->filenameMUT,  "w");
  metil_file->STAT = fopen(metil_file->filenameSTAT, "w");


  FILE *a;
  a = metil_file->CpG;
  file_error = fprintf(a, "File for Cytosines in CpG context\n\n");
  if (file_error < 0) {
    printf("Error al escribir\n");
    exit(-1);
  }
  a = metil_file->CHG;
  file_error = fprintf(a, "File for Cytosines in CHG context\n\n");
  if (file_error < 0) {
    printf("Error al escribir\n");
    exit(-1);
  }
  a = metil_file->CHH;
  file_error = fprintf(a, "File for Cytosines in CHH context\n\n");
  if (file_error < 0) {
    printf("Error al escribir\n");
    exit(-1);
  }
  a = metil_file->MUT;
  file_error = fprintf(a, "File for Cytosines mutated\n\n");
  if (file_error < 0) {
    printf("Error al escribir\n");
    exit(-1);
  }
  a = metil_file->STAT;
  file_error = fprintf(a, "File for Methylation Statistics\n\n");
  if (file_error < 0) {
    printf("Error al escribir\n");
    exit(-1);
  }

  free(name_tmp);

  metil_file->CpG_methyl   = 0;
  metil_file->CpG_unmethyl = 0;
  metil_file->CHG_methyl   = 0;
  metil_file->CHG_unmethyl = 0;
  metil_file->CHH_methyl   = 0;
  metil_file->CHH_unmethyl = 0;
  metil_file->MUT_methyl   = 0;
  metil_file->num_bases    = 0;

  /*
  metil_file->table_isles = hash_create(jenkins_one_at_a_time_hash,
					scmp,
					destr_key,
					nulldes,
					200);
  hash_init(metil_file->table_isles, "ACGTN-");
  */
}

//====================================================================================

void metil_file_free(metil_file_t *metil_file) {
  free(metil_file->filenameCpG);
  free(metil_file->filenameCHG);
  free(metil_file->filenameCHH);
  free(metil_file->filenameMUT);
  free(metil_file->filenameSTAT);

  if (metil_file->CpG  != NULL) fclose(metil_file->CpG);
  if (metil_file->CHG  != NULL) fclose(metil_file->CHG);
  if (metil_file->CHH  != NULL) fclose(metil_file->CHH);
  if (metil_file->MUT  != NULL) fclose(metil_file->MUT);
  if (metil_file->STAT != NULL) fclose(metil_file->STAT);

  /*
  hash_destroy(metil_file->table_isles);
  */
  free(metil_file);
}

//====================================================================================

void write_metilation_status(array_list_t *array_list, metil_file_t *metil_file) {

  //printf("Init metilation status\n");
  
  FILE * CpG = metil_file->CpG;
  if (CpG == NULL) {
    printf("reopen CpG file\n");
    CpG =fopen(metil_file->filenameCpG, "a");
  }
  FILE * CHG = metil_file->CHG;
  if (CHG == NULL) {
    printf("reopen CHG file\n");
    CHG = fopen(metil_file->filenameCHG, "a");
  }
  FILE * CHH = metil_file->CHH;
  if (CHH == NULL) {
    printf("reopen CHH file\n");
    CHH = fopen(metil_file->filenameCHH, "a");
  }
  FILE * MUT = metil_file->MUT;
  if (MUT == NULL) {
    printf("reopen CHH file\n");
    MUT = fopen(metil_file->filenameMUT, "a");
  }

  size_t num_items = array_list_size(array_list);
  alignment_t *alig;
  genome_t *genome = metil_file->genome;
  hash_table_t *table_bs = metil_file->table_isles;
  char *seq, *gen;
  size_t len, end, start;
  char *key = (char *)malloc(3 * sizeof(char));
  int data;

  size_t contador = 0;
  size_t alineamientos = 0;

  for (size_t j = 0; j < num_items; j++) {
    alig = (alignment_t *) array_list_get(j, array_list);
    if (alig != NULL && alig->is_seq_mapped) {
      /*
      printf("alignment %lu\n", j);
      printf("query_name %s\n", alig->query_name);
      printf("sequence   %s\n", alig->sequence);
      printf("cigar      %s\n", alig->cigar);
      printf("position   %i\n", alig->position);
      printf("mate pos   %i\n", alig->mate_position);
      printf("temp len   %i\n", alig->template_length);
      printf("chromo     %i\n", alig->chromosome);
      printf("strand     %i\n", alig->seq_strand);
      printf("mate  chro %i\n", alig->mate_chromosome);
      printf("map qual   %i\n", alig->map_quality);
      printf("n cigar    %i\n", alig->num_cigar_operations);
      */      
      
      //seq = obtain_seq(alig);
      //printf("seq %s\n", seq);

      len = strlen(seq);
      gen = (char *)calloc(len + 2, sizeof(char));

      start = alig->position;
      end = start + len - 1;
      if (end >= genome->chr_size[alig->chromosome]) {
	//printf("    end %lu\n", end);
	//printf("    gen %lu\n", genome->chr_size[alig->chromosome]);
        end = genome->chr_size[alig->chromosome] - 1;
	//printf("new end %lu\n", end);
      }


      //printf("seq %s\n", seq);
      //printf("chromo %i, strand %i, begin %lu, end %lu\n",
      //     alig->chromosome, alig->seq_strand, alig->position, end);
      
      genome_read_sequence_by_chr_index(gen, alig->seq_strand, alig->chromosome, &start, &end, genome);

      //printf("\nseq %s\ngen %s\n", seq, gen);
      
      alineamientos++;
      for (size_t i = 0; i < len - 2; i++) {
	if (gen[i] == 'C' && seq[i] == 'C') contador++;
	if (gen[i] == 'G' && seq[i] == 'G') contador++;
	/*
	if (gen[i] == 'C') {
	  if (seq[i] == 'C' || seq[i] == 'T') {
	    key[0] = gen[i];
	    key[1] = gen[i + 1];
	    key[2] = gen[i + 2];
	    data = (int)hash_find_data(table_bs, key);
	    
	    //printf("Candidata (%i) en %lu\n", data, start + i);

	    switch (data) {
	    case ZONE_CpG:
	      if (seq[i] == 'C') {
		printf("%s\t+\t%i %i\t%lu\tZ\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      } else {
		printf("%s\t-\t%i %i\t%lu\tz\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      }
	      break;
	    case ZONE_CHG:
	      if (seq[i] == 'C') {
		printf("%s\t+\t%i %i\t%lu\tX\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      } else {
		printf("%s\t-\t%i %i\t%lu\tx\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      }
	      break;
	    case ZONE_CHH:
	      if (seq[i] == 'C') {
		printf("%s\t+\t%i %i\t%lu\tH\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      } else {
		printf("%s\t-\t%i %i\t%lu\th\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      }
	      break;
	    case ZONE_OTHER:
	      if (seq[i] == 'C') {
		printf("%s\t+\t%i %i\t%lu\tM\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      } else {
		printf("%s\t-\t%i %i\t%lu\tm\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      }
	      break;
	    default:
	      printf("%s\t-\t%i %i\t%lu\t???\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	    }
	  }
	  else {
	    printf("%s\t-\t%i %i\t%lu\tm\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	    //printf("Mutacion en %lu\n", start + i);
	  }
	}
	*/
      }

      if (seq) free(seq);
      if (gen) free(gen);
    }
  }
  
  printf("methylated cytosines:\t%lu\tin\t%lu\talignments\n", contador, alineamientos);
  
  free(key);
}

//====================================================================================

char *obtain_seq_old(alignment_t *alig) {
  char *read = alig->sequence;
  char *cigar = strdup(alig->cigar);
  //char *cigar = strdup("100M1D");
  int num;
  char car;
  int cont, pos, pos_read;
  int len = strlen(cigar) - 1;
  char *seq = (char *)calloc(1024, sizeof(char));
  //printf("cigar %s\n", cigar);
  //printf("read  %s\n", read);
  
  pos = 0;
  pos_read = 0;
  while(strlen(cigar) > 0 && strlen(cigar) != len) {
    len = strlen(cigar);
    sscanf(cigar, "%i%c%s", &num, &car, cigar);
    //printf("%3i %c %s\n",  num, car, cigar);
    if (car == 'M' || car == '=') {
      for (cont = 0; cont < num; cont++, pos++, pos_read++) {
	seq[pos] = read[pos_read];
      }
    }
    else {
      if (car == 'D') {
	pos_read += num - 1;
      }
      else {
	if (car == 'I') {
	  for (cont = 0; cont < num; cont++, pos++) {
	    seq[pos] = 'N';
	  }
	}
      }
    }
  }
  seq[pos] = '\0';

  //printf("seq   %s\n", seq);

  free(cigar);
  return seq;
}

//====================================================================================

void remove_duplicates(size_t reads, array_list_t **list, array_list_t **list2) {
  size_t num_items, num_items2;
  alignment_t *alig, *alig2;

  for (size_t i = 0; i < reads; i++) {
    num_items = array_list_size(list[i]);
    //printf("list[%lu]\tlist2[%lu]\n", array_list_size(list[i]), array_list_size(list2[i]));
    for (size_t j = 0; j < num_items; j++) {
      //printf("list[%lu]\tlist2[%lu]\n", array_list_size(list[i]), array_list_size(list2[i]));
      alig = (alignment_t *) array_list_get(j, list[i]);

      if (alig != NULL && alig->is_seq_mapped) {
	num_items2 = array_list_size(list2[i]);
	for (size_t k = 0; k < num_items2; k++) {
	//for (size_t k = num_items2 - 1; k >= 0; k--) {
	  alig2 = (alignment_t *) array_list_get(k, list2[i]);
	  /*
	  printf("alignment %lu - %lu\n", j, k);
	  printf("query_name %s\n",  alig->query_name);
	  printf("query_name %s\n", alig2->query_name);
	  printf("sequence   %s\n",  alig->sequence);
	  printf("sequence   %s\n", alig2->sequence);
	  printf("cigar      %s\n",  alig->cigar);
	  printf("cigar      %s\n", alig2->cigar);
	  printf("position   %i\n",  alig->position + 1);
	  printf("position   %i\n", alig2->position + 1);
	  printf("mate pos   %i\n",  alig->mate_position);
	  printf("mate pos   %i\n", alig2->mate_position);
	  printf("temp len   %i\n",  alig->template_length);
	  printf("temp len   %i\n", alig2->template_length);
	  printf("chromo     %i\n",  alig->chromosome);
	  printf("chromo     %i\n", alig2->chromosome);
	  printf("strand     %i\n",  alig->seq_strand);
	  printf("strand     %i\n", alig2->seq_strand);
	  printf("mate  chro %i\n",  alig->mate_chromosome);
	  printf("mate  chro %i\n", alig2->mate_chromosome);
	  printf("map qual   %i\n",  alig->map_quality);
	  printf("map qual   %i\n", alig2->map_quality);
	  printf("n cigar    %i\n",  alig->num_cigar_operations);
	  printf("n cigar    %i\n", alig2->num_cigar_operations);
	  */
	  if (alig->position   == alig2->position
	      && alig->chromosome == alig2->chromosome 
	      && alig->seq_strand == alig2->seq_strand 
	      && alig->num_cigar_operations == alig2->num_cigar_operations) {
	    //printf("list[%lu]\tlist2[%lu]\n", array_list_size(list[i]), array_list_size(list2[i]));
	    alig2 = (alignment_t *) array_list_remove_at(k, list2[i]);
	    alignment_free(alig2);
	    //printf("list[%lu]\tlist2[%lu]\n", array_list_size(list[i]), array_list_size(list2[i]));
	    k--;
	    num_items2 = array_list_size(list2[i]);
	  }
	}
      }
    }
    //printf("list[%lu]\tlist2[%lu]\n", array_list_size(list[i]), array_list_size(list2[i]));
  }
}

//====================================================================================

int methylation_status_report(sw_server_input_t* input, batch_t *batch) {

  mapping_batch_t *mapping_batch = (mapping_batch_t *) batch->mapping_batch;
  array_list_t **mapping_lists;
  size_t num_items;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  genome_t *genome = input->genome_p;
  fastq_read_t *orig;
  

  // inicializar listas para guardar datos de c's metiladas/no metiladas
  bs_context_t *bs_context = bs_context_new(10000);
  mapping_batch->bs_context = bs_context;

  /*
  printf("Lists:\nCpG %lu\tCHG %lu\tCHH %lu\tCMUT %lu\n",
	 array_list_size(bs_context->context_CpG), array_list_size(bs_context->context_CHG),
	 array_list_size(bs_context->context_CHH), array_list_size(bs_context->context_MUT));
  printf("\tCpG\tCHG\tCHH\tMUT\nMet\t%lu\t%lu\t%lu\t%lu\nUnMet\t%lu\t%lu\t%lu\nBases\t%lu\n",
	 bs_context->CpG_methyl,
	 bs_context->CHG_methyl,
	 bs_context->CHH_methyl,
	 bs_context->MUT_methyl,
	 bs_context->CpG_unmethyl,
	 bs_context->CHG_unmethyl,
	 bs_context->CHH_unmethyl,
	 bs_context->num_bases);
  */

  remove_duplicates(num_reads, mapping_batch->mapping_lists, mapping_batch->mapping_lists2);
  
  for (int k = 0; k < 2; k++) {
    mapping_lists = (k == 0) ? mapping_batch->mapping_lists : mapping_batch->mapping_lists2;
    
    for (size_t i = 0; i < num_reads; i++) {
      num_items = array_list_size(mapping_lists[i]);
      
      // mapped or not mapped ?
      if (num_items != 0) {
	add_metilation_status(mapping_lists[i], bs_context, genome, mapping_batch->fq_batch, i, k);
	//add_metilation_status_bin(mapping_lists[i], bs_context, input->valuesCT, input->valuesGA, mapping_batch->fq_batch, i, k);
      }
    }
  }

  //  if (array_list_size(bs_context->context_CpG) > 10000)
  /*
  printf("----------------> size = %lu %lu %lu %lu\n", 
	 array_list_size(bs__context->context_CpG),
	 array_list_size(bs_context->context_CpG), 
	 array_list_size(bs_context->context_CpG), 
	 array_list_size(bs_context->context_MUT));
  */
  return CONSUMER_STAGE;
}

//====================================================================================

void add_metilation_status(array_list_t *array_list, bs_context_t *bs_context, genome_t * genome, array_list_t * orig_seq, size_t index, int conversion) {

  //printf("Init add metilation status\n");

  size_t num_items = array_list_size(array_list);
  alignment_t *alig;
  char *seq, *gen;
  fastq_read_t *orig;
  size_t len, end, start;
  int new_strand;
  char *new_stage;
  metil_data_t *metil_data;

  int write_file = 1;

  orig = (fastq_read_t *) array_list_get(index, orig_seq);

  for (size_t j = 0; j < num_items; j++) {

    alig = (alignment_t *) array_list_get(j, array_list);

    if (alig != NULL && alig->is_seq_mapped) {

      /*
      printf("alignment %lu\n", j);
      printf("query_name %s\n", alig->query_name);
      printf("sequence   %s\n", alig->sequence);
      printf("cigar      %s\n", alig->cigar);
      printf("position   %i\n", alig->position + 1);
      printf("mate pos   %i\n", alig->mate_position);
      printf("temp len   %i\n", alig->template_length);
      printf("chromo     %i\n", alig->chromosome);
      printf("strand     %i\n", alig->seq_strand);
      printf("mate  chro %i\n", alig->mate_chromosome);
      printf("map qual   %i\n", alig->map_quality);
      printf("n cigar    %i\n", alig->num_cigar_operations);
      */

  /**************
  comprobar esta funcion para optimizarla
  **************/
      seq = obtain_seq(alig, orig);

      if (alig->seq_strand == 1) {
	char *seq_dup = strdup(seq);
	rev_comp(seq_dup, seq, orig->length);
	free(seq_dup);
      }

      // increase de counter number of bases
      //bs_context->num_bases += orig->length;

      len = orig->length;
      gen = (char *)calloc(len + 5, sizeof(char));

      start = alig->position + 1;
      end = start + len + 4;
      if (end >= genome->chr_size[alig->chromosome]) {
        end = genome->chr_size[alig->chromosome] - 1;
      }

      genome_read_sequence_by_chr_index(gen, alig->seq_strand, alig->chromosome, &start, &end, genome);

      /*
      printf("seq %s\n", seq);
      printf("gen %s\n", gen);
      */
      /*
      printf("chromo %i, strand %i, begin %lu, end %lu\n",
	     alig->chromosome, alig->seq_strand, alig->position, end);
      */


      for (size_t i = 0; i < len; i++) {
	if ((conversion == 1 && alig->seq_strand == 0) || (conversion == 0 && alig->seq_strand == 1)) {
	  // methylated/unmethylated cytosines are located in the same strand as the alignment
	  if (gen[i] == 'C') {
	    //case ZONE_CpG:
	    if (gen[i + 1] == 'G') {
	      if (seq[i] == 'C') {
		if (write_file == 1)
		  postprocess_bs(alig->query_name, '+', alig->chromosome, start + i, 'Z', alig->seq_strand, 0,
				 bs_context->context_CpG);
		bs_context->CpG_methyl++;
	      } else if (seq[i] == 'T') {
		if (write_file == 1)
		  postprocess_bs(alig->query_name, '-', alig->chromosome, start + i, 'z', alig->seq_strand, 0,
				 bs_context->context_CpG);
		bs_context->CpG_unmethyl++;
	      } else {
		if (write_file == 1)
		  postprocess_bs(alig->query_name, '.', alig->chromosome, start + i, 'M', alig->seq_strand, 3,
				 bs_context->context_MUT);
		bs_context->MUT_methyl++;
	      }
	    } else {
	      //case ZONE_CHG:
	      if (gen[i + 2] == 'G') {
		if (seq[i] == 'C') {
		  if (write_file == 1)
		    postprocess_bs(alig->query_name, '+', alig->chromosome, start + i, 'X', alig->seq_strand, 1,
				   bs_context->context_CHG);
		  bs_context->CHG_methyl++;
		} else if (seq[i] == 'T') {
		  if (write_file == 1)
		    postprocess_bs(alig->query_name, '-', alig->chromosome, start + i, 'x', alig->seq_strand, 1,
				   bs_context->context_CHG);
		  bs_context->CHG_unmethyl++;
		} else {
		  if (write_file == 1)
		    postprocess_bs(alig->query_name, '.', alig->chromosome, start + i, 'M', alig->seq_strand, 3,
				   bs_context->context_MUT);
		  bs_context->MUT_methyl++;
		}
	      } else {
		//case ZONE_CHH:
		if (seq[i] == 'C') {
		  if (write_file == 1)
		    postprocess_bs(alig->query_name, '+', alig->chromosome, start + i, 'H', alig->seq_strand, 2,
				   bs_context->context_CHH);
		  bs_context->CHH_methyl++;
		} else if (seq[i] == 'T') {
		  if (write_file == 1)
		    postprocess_bs(alig->query_name, '-', alig->chromosome, start + i, 'h', alig->seq_strand, 2,
				   bs_context->context_CHH);
		  bs_context->CHH_unmethyl++;
		} else {
		  if (write_file == 1)
		    postprocess_bs(alig->query_name, '.', alig->chromosome, start + i, 'M', alig->seq_strand, 3,
				   bs_context->context_MUT);
		  bs_context->MUT_methyl++;
		}
	      }
	    }
	  }
	} else {
	  // methylated/unmethylated cytosines are located in the other strand
	  if (alig->seq_strand == 0) new_strand = 1;
	  else                       new_strand = 0;
	  
	  if (gen[i+2] == 'G') {
	    //case ZONE_CpG:
	    if (gen[i + 1] == 'C') {
	      if (seq[i] == 'G') {
		if (write_file == 1)
		  postprocess_bs(alig->query_name, '+', alig->chromosome, start + i, 'Z', new_strand, 0,
				 bs_context->context_CpG);
		bs_context->CpG_methyl++;
	      } else if (seq[i] == 'A') {
		if (write_file == 1)
		  postprocess_bs(alig->query_name, '-', alig->chromosome, start + i, 'z', new_strand, 0,
			       bs_context->context_CpG);
		bs_context->CpG_unmethyl++;
	      } else {
		if (write_file == 1)
		  postprocess_bs(alig->query_name, '.', alig->chromosome, start + i, 'M', alig->seq_strand, 3,
				 bs_context->context_MUT);
		bs_context->MUT_methyl++;
	      }
	    } else {
	      //case ZONE_CHG:
	      if (gen[i] == 'C') {
		if (seq[i] == 'G') {
		  if (write_file == 1)
		    postprocess_bs(alig->query_name, '+', alig->chromosome, start + i, 'X', new_strand, 1,
				   bs_context->context_CHG);
		  bs_context->CHG_methyl++;
		} else if (seq[i] == 'A') {
		  if (write_file == 1)
		    postprocess_bs(alig->query_name, '-', alig->chromosome, start + i, 'x', new_strand, 1,
				   bs_context->context_CHG);
		  bs_context->CHG_unmethyl++;
		} else {
		  if (write_file == 1)
		    postprocess_bs(alig->query_name, '.', alig->chromosome, start + i, 'M', alig->seq_strand, 3,
				   bs_context->context_MUT);
		  bs_context->MUT_methyl++;
		}
	      } else {
		//case ZONE_CHH:
		if (seq[i] == 'G') {
		  if (write_file == 1)
		    postprocess_bs(alig->query_name, '+', alig->chromosome, start + i, 'H', new_strand, 2,
				   bs_context->context_CHH);
		  bs_context->CHH_methyl++;
		} else if (seq[i] == 'A') {
		  if (write_file == 1)
		    postprocess_bs(alig->query_name, '-', alig->chromosome, start + i, 'h', new_strand, 2,
				   bs_context->context_CHH);
		  bs_context->CHH_unmethyl++;
		} else {
		  if (write_file == 1)
		    postprocess_bs(alig->query_name, '.', alig->chromosome, start + i, 'M', alig->seq_strand, 3,
				   bs_context->context_MUT);
		  bs_context->MUT_methyl++;
		}
	      }
	    }
	  }
	}
	//printf("End for\n");
      }

      if (seq) free(seq);
      if (gen) free(gen);

    }

  }

  /*
  printf("\tMethyl\tunMethyl\nCpG\t%lu\t%lu\nCHG\t%lu\t%lu\nCHH\t%lu\t%lu\n------------------------\nMUT\t%lu\n",
	 bs_context->CpG_methyl, bs_context->CpG_unmethyl,
	 bs_context->CHG_methyl, bs_context->CHG_unmethyl,
	 bs_context->CHH_methyl, bs_context->CHH_unmethyl,
	 bs_context->MUT_methyl);
  */
}

//====================================================================================

void add_metilation_status_bin(array_list_t *array_list, bs_context_t *bs_context,
			       unsigned long long **gen_binCT, unsigned long long **gen_binGA,
			       array_list_t * orig_seq, size_t index, int conversion) {

  //printf("Init add metilation status\n");

  size_t num_items = array_list_size(array_list);
  alignment_t *alig;
  char *seq, *gen;
  fastq_read_t *orig;
  size_t len, end, start;
  int new_strand;
  char *new_stage;
  metil_data_t *metil_data;
  char posit = '+', negat = '-';

  int write_file = 1;

  orig = (fastq_read_t *) array_list_get(index, orig_seq);

  for (size_t j = 0; j < num_items; j++) {

    alig = (alignment_t *) array_list_get(j, array_list);

    if (alig != NULL && alig->is_seq_mapped) {

  /**************
  comprobar esta funcion para optimizarla
  **************/
      seq = obtain_seq(alig, orig);

      if (alig->seq_strand == 1) {
	char *seq_dup = strdup(seq);
	rev_comp(seq_dup, seq, orig->length);
	free(seq_dup);
      }

      len = orig->length;

      start = alig->position;
      end = start + len;
/*
      if (end >= genome->chr_size[alig->chromosome]) {
        end = genome->chr_size[alig->chromosome] - 1;
      }
*/
      //printf("++++++++ CT = %llu -------- GA = %llu\n", gen_binCT[0][100000], gen_binGA[0][100000]);

      if ((conversion == 1 && alig->seq_strand == 0) || (conversion == 0 && alig->seq_strand == 1)) {
        search_methylation(alig->chromosome, start, end, gen_binCT, seq, bs_context, '+', 0, alig->query_name);
      }
      else {
        search_methylation(alig->chromosome, start, end, gen_binCT, seq, bs_context, '-', 1, alig->query_name);
      }

      if (seq) free(seq);
      if (gen) free(gen);

    }
  }

  /*
  printf("\tMethyl\tunMethyl\nCpG\t%lu\t%lu\nCHG\t%lu\t%lu\nCHH\t%lu\t%lu\n------------------------\nMUT\t%lu\n",
	 bs_context->CpG_methyl, bs_context->CpG_unmethyl,
	 bs_context->CHG_methyl, bs_context->CHG_unmethyl,
	 bs_context->CHH_methyl, bs_context->CHH_unmethyl,
	 bs_context->MUT_methyl);
  */
}

//====================================================================================

void write_bs_context(metil_file_t *metil_file, bs_context_t *bs_context) {

  array_list_t *context_CpG = bs_context->context_CpG;
  array_list_t *context_CHG = bs_context->context_CHG;
  array_list_t *context_CHH = bs_context->context_CHH;
  array_list_t *context_MUT = bs_context->context_MUT;

  size_t num_items, num_reads;
  char *bs_seq;
  int file_error;
  metil_data_t *metil_data;

  FILE * CpG = metil_file->CpG;
  FILE * CHG = metil_file->CHG;
  FILE * CHH = metil_file->CHH;
  FILE * MUT = metil_file->MUT;

  /*
  FILE * CpG;
  FILE * CHG;
  FILE * CHH;
  FILE * MUT;
  */

  metil_file->CpG_methyl   += bs_context->CpG_methyl;
  metil_file->CpG_unmethyl += bs_context->CpG_unmethyl;
  metil_file->CHG_methyl   += bs_context->CHG_methyl;
  metil_file->CHG_unmethyl += bs_context->CHG_unmethyl;
  metil_file->CHH_methyl   += bs_context->CHH_methyl;
  metil_file->CHH_unmethyl += bs_context->CHH_unmethyl;
  metil_file->MUT_methyl   += bs_context->MUT_methyl;
  metil_file->num_bases    += bs_context->num_bases;

  if (CpG == NULL) {
    printf("reopen CpG file\n");
    CpG =fopen(metil_file->filenameCpG, "a");
  }
  if (CHG == NULL) {
    printf("reopen CHG file\n");
    CHG = fopen(metil_file->filenameCHG, "a");
  }
  if (CHH == NULL) {
    printf("reopen CHH file\n");
    CHH = fopen(metil_file->filenameCHH, "a");
  }
  if (MUT == NULL) {
    printf("reopen CHH file\n");
    MUT = fopen(metil_file->filenameMUT, "a");
  }

  if (context_CpG != NULL) {
    num_items = array_list_size(context_CpG);
    for (int i = num_items - 1; i >= 0; i--) {
      //for (size_t i = 0; i < num_items; i++) {
      metil_data = (metil_data_t *)array_list_get(i, context_CpG);
      file_error = fprintf(CpG, "%s\t%c\t%i\t%lu\t%c\t%i\n",
			   metil_data->query_name, metil_data->status,
			   metil_data->chromosome, metil_data->start,
			   metil_data->context, metil_data->strand);
      metil_data_free(metil_data);
      
      if (file_error < 0) {
	printf("Error on write\n");
	exit(-1);
      }
    }
  }
  //if (context_CpG) array_list_free(context_CpG, NULL);
  
  if (context_CHG != NULL) {
    num_items = array_list_size(context_CHG);
    for (int i = num_items - 1; i >= 0; i--) {
      //for (size_t i = 0; i < num_items; i++) {
      metil_data = (metil_data_t *)array_list_get(i, context_CHG);
      file_error = fprintf(CHG, "%s\t%c\t%i\t%lu\t%c\t%i\n",
			   metil_data->query_name, metil_data->status,
			   metil_data->chromosome, metil_data->start,
			   metil_data->context, metil_data->strand);
      metil_data_free(metil_data);
      
      if (file_error < 0) {
	printf("Error on write\n");
	exit(-1);
      }
    }
  }
  //if (context_CHG) array_list_free(context_CHG, NULL);
  
  if (context_CHH != NULL) {
    num_items = array_list_size(context_CHH);
    for (int i = num_items - 1; i >= 0; i--) {
      //for (size_t i = 0; i < num_items; i++) {
      metil_data = (metil_data_t *)array_list_get(i, context_CHH);
      file_error = fprintf(CHH, "%s\t%c\t%i\t%lu\t%c\t%i\n",
			   metil_data->query_name, metil_data->status,
			   metil_data->chromosome, metil_data->start,
			   metil_data->context, metil_data->strand);
      metil_data_free(metil_data);
      
      if (file_error < 0) {
	printf("Error on write\n");
	exit(-1);
      }
    }
  }
  //if (context_CHH) array_list_free(context_CHH, NULL);
  
  if (context_MUT != NULL) {
    num_items = array_list_size(context_MUT);
    for (int i = num_items - 1; i >= 0; i--) {
      //for (size_t i = 0; i < num_items; i++) {
      metil_data = (metil_data_t *)array_list_get(i, context_MUT);
      file_error = fprintf(MUT, "%s\t%c\t%i\t%lu\t%c\t%i\n",
			   metil_data->query_name, metil_data->status,
			   metil_data->chromosome, metil_data->start,
			   metil_data->context, metil_data->strand);
      metil_data_free(metil_data);
      
      if (file_error < 0) {
	printf("Error on write\n");
	exit(-1);
      }
    }
  }
  //if (context_MUT) array_list_free(context_MUT, NULL);
  
  if (context_CpG) array_list_free(context_CpG, NULL);
  if (context_CHG) array_list_free(context_CHG, NULL);
  if (context_CHH) array_list_free(context_CHH, NULL);
  if (context_MUT) array_list_free(context_MUT, NULL);
  if (bs_context)  bs_context_free(bs_context);

  /*  
  // prueba de escritura
  fclose(CpG);
  // fin prueba de escritura
  */
  /*
  printf("\nWriter\n\tMethyl\tunMethyl\nCpG\t%lu\t%lu\nCHG\t%lu\t%lu\nCHH\t%lu\t%lu\n------------------------\nMUT\t%lu\n",
	 metil_file->CpG_methyl, metil_file->CpG_unmethyl,
	 metil_file->CHG_methyl, metil_file->CHG_unmethyl,
	 metil_file->CHH_methyl, metil_file->CHH_unmethyl,
	 metil_file->MUT_methyl);
  */
  /*
  metil_file->num_bases += 1;
  printf("\nAligner\n\tMethyl\tunMethyl\nCpG\t%7.2f\t%7.2f\nCHG\t%7.2f\t%7.2f\nCHH\t%7.2f\t%7.2f\n------------------------\nMUT\t%7.2\
f\n",
         (float) metil_file->CpG_methyl / metil_file->num_bases * 100,
         (float) metil_file->CpG_unmethyl / metil_file->num_bases * 100,
         (float) metil_file->CHG_methyl / metil_file->num_bases * 100,
         (float) metil_file->CHG_unmethyl / metil_file->num_bases * 100,
         (float) metil_file->CHH_methyl / metil_file->num_bases * 100,
         (float) metil_file->CHH_unmethyl / metil_file->num_bases * 100,
         (float) metil_file->MUT_methyl / metil_file->num_bases * 100);
  */
}

//====================================================================================

void write_context_bs(metil_file_t *metil_file, bs_context_t *bs_context) {

  array_list_bs_t *context_CpG = bs_context->context_bs_CpG;
  array_list_bs_t *context_CHG = bs_context->context_bs_CHG;
  array_list_bs_t *context_CHH = bs_context->context_bs_CHH;
  array_list_bs_t *context_MUT = bs_context->context_bs_MUT;

  size_t num_items, num_reads;
  char *bs_seq;
  int file_error;
  metil_data_t *metil_data;

  FILE * CpG = metil_file->CpG;
  FILE * CHG = metil_file->CHG;
  FILE * CHH = metil_file->CHH;
  FILE * MUT = metil_file->MUT;

  metil_file->CpG_methyl   += bs_context->CpG_methyl;
  metil_file->CpG_unmethyl += bs_context->CpG_unmethyl;
  metil_file->CHG_methyl   += bs_context->CHG_methyl;
  metil_file->CHG_unmethyl += bs_context->CHG_unmethyl;
  metil_file->CHH_methyl   += bs_context->CHH_methyl;
  metil_file->CHH_unmethyl += bs_context->CHH_unmethyl;
  metil_file->MUT_methyl   += bs_context->MUT_methyl;
  metil_file->num_bases    += bs_context->num_bases;

  if (CpG == NULL) {
    printf("reopen CpG file\n");
    CpG =fopen(metil_file->filenameCpG, "a");
  }
  if (CHG == NULL) {
    printf("reopen CHG file\n");
    CHG = fopen(metil_file->filenameCHG, "a");
  }
  if (CHH == NULL) {
    printf("reopen CHH file\n");
    CHH = fopen(metil_file->filenameCHH, "a");
  }
  if (MUT == NULL) {
    printf("reopen CHH file\n");
    MUT = fopen(metil_file->filenameMUT, "a");
  }
  
  if (context_CpG != NULL) {
    num_items = array_list_bs_size(context_CpG);
    for (int i = num_items - 1; i >= 0; i--) {
      //for (size_t i = 0; i < num_items; i++) {
      metil_data = (metil_data_t *)array_list_bs_get(i, context_CpG);
      file_error = fprintf(CpG, "%s\t%c\t%i\t%lu\t%c\t%i\n",
			   metil_data->query_name, metil_data->status,
			   metil_data->chromosome, metil_data->start,
			   metil_data->context, metil_data->strand);
      //metil_data_free(metil_data);
      
      if (file_error < 0) {
	printf("Error on write\n");
	exit(-1);
      }
    }
  }
  //if (context_CpG) array_list_free(context_CpG, NULL);
  
  if (context_CHG != NULL) {
    num_items = array_list_bs_size(context_CHG);
    for (int i = num_items - 1; i >= 0; i--) {
      //for (size_t i = 0; i < num_items; i++) {
      metil_data = (metil_data_t *)array_list_bs_get(i, context_CHG);
      file_error = fprintf(CHG, "%s\t%c\t%i\t%lu\t%c\t%i\n",
			   metil_data->query_name, metil_data->status,
			   metil_data->chromosome, metil_data->start,
			   metil_data->context, metil_data->strand);
      //metil_data_free(metil_data);
      
      if (file_error < 0) {
	printf("Error on write\n");
	exit(-1);
      }
    }
  }
  //if (context_CHG) array_list_free(context_CHG, NULL);
  
  if (context_CHH != NULL) {
    num_items = array_list_bs_size(context_CHH);
    for (int i = num_items - 1; i >= 0; i--) {
      //for (size_t i = 0; i < num_items; i++) {
      metil_data = (metil_data_t *)array_list_bs_get(i, context_CHH);
      file_error = fprintf(CHH, "%s\t%c\t%i\t%lu\t%c\t%i\n",
			   metil_data->query_name, metil_data->status,
			   metil_data->chromosome, metil_data->start,
			   metil_data->context, metil_data->strand);
      //metil_data_free(metil_data);
      
      if (file_error < 0) {
	printf("Error on write\n");
	exit(-1);
      }
    }
  }
  //if (context_CHH) array_list_free(context_CHH, NULL);
  
  if (context_MUT != NULL) {
    num_items = array_list_bs_size(context_MUT);
    for (int i = num_items - 1; i >= 0; i--) {
      //for (size_t i = 0; i < num_items; i++) {
      metil_data = (metil_data_t *)array_list_bs_get(i, context_MUT);
      file_error = fprintf(MUT, "%s\t%c\t%i\t%lu\t%c\t%i\n",
			   metil_data->query_name, metil_data->status,
			   metil_data->chromosome, metil_data->start,
			   metil_data->context, metil_data->strand);
      //metil_data_free(metil_data);
      
      if (file_error < 0) {
	printf("Error on write\n");
	exit(-1);
      }
    }
  }
  //if (context_MUT) array_list_free(context_MUT, NULL);

  if (context_CpG) array_list_bs_free(context_CpG);
  if (context_CHG) array_list_bs_free(context_CHG);
  if (context_CHH) array_list_bs_free(context_CHH);
  if (context_MUT) array_list_bs_free(context_MUT);
  if (bs_context)  bs_context_free(bs_context);

  /*
  printf("\nWriter\n\tMethyl\tunMethyl\nCpG\t%lu\t%lu\nCHG\t%lu\t%lu\nCHH\t%lu\t%lu\n------------------------\nMUT\t%lu\n",
	 metil_file->CpG_methyl, metil_file->CpG_unmethyl,
	 metil_file->CHG_methyl, metil_file->CHG_unmethyl,
	 metil_file->CHH_methyl, metil_file->CHH_unmethyl,
	 metil_file->MUT_methyl);
  */
}

//====================================================================================

void metil_data_init(metil_data_t *metil_data, char *query, char status, int chromosome, size_t start, char context, int strand, int zone) {
  if (metil_data == NULL)
    metil_data = (metil_data_t *)malloc(sizeof(metil_data_t));

  metil_data->query_name = strdup(query);
  metil_data->status = status;
  metil_data->chromosome = chromosome;
  metil_data->start = start;
  metil_data->context = context;
  metil_data->strand = strand;
  metil_data->zone = zone;
}

//====================================================================================

void metil_data_free(metil_data_t *metil_data) {
  if (metil_data != NULL) {
    free(metil_data->query_name);
    free(metil_data);
  }
}

//====================================================================================

void postprocess_bs(char *query_name, char status, size_t chromosome, size_t start, char context, int strand, int region,
		    array_list_t *list) {

  metil_data_t *metil_data = (metil_data_t *) malloc(sizeof(metil_data_t));
  metil_data_init(metil_data, query_name, status, chromosome, start, context, strand, region);
  array_list_insert(metil_data, list);

}

//====================================================================================

void postproc_bs(char *query_name, char status, size_t chromosome, size_t start, char context, int strand, int region,
		 array_list_bs_t *list) {
  array_list_bs_insert(list, query_name, status, chromosome, start, context, strand, region);
  //printf("Entro\n");
}

//====================================================================================

void search_methylation(int c, size_t init, size_t end, unsigned long long **values, char *read, bs_context_t *bs_context, char strand, int type, char *query_name){
  size_t init_num, end_num;
  int i, init_num_pos, end_num_pos;
  unsigned long long tmp, res;
  size_t j, pos;
  size_t algo = 0;

  int write_file = 1;

  int elem = (sizeof(unsigned long long) << 2);
  //printf("elem = %i\n", elem);

  init_num = init / elem;
  init_num_pos = init % elem;
  end_num = end / elem;
  end_num_pos = end % elem;
  char c1, c2;
  if (type == 0) {
    c1 = 'C'; c2 = 'T';
  } else {
    c1 = 'G'; c2 = 'A';
  }
  
  /*
  printf("Chromosome %2i (%10lu - %10lu)\n", c, init, end);
  printf("From number\t= %10lu (%2i)\nTo number\t= %10lu (%2i)\n\n",
	 init_num, init_num_pos, end_num, end_num_pos);
  printf("value = %llu\n", values[c][init_num]);
  */

  // recorrer los valores de derecha a izquierda
  tmp  = values[c][end_num];

  // descartar los pares de bits menos significativos hasta llegar al valor donde termina la cadena
  for (i = elem; i > end_num_pos; i--) {
    tmp = tmp >> 2;
    //tmp /= 4;
  }


  // recorrer los pares de bits comprobando el contexto del ultimo numero
  for (i = end_num_pos, pos = 0; i >= 0; i--, pos++) {
    search_gen();
    //search_gen(tmp, read, pos, c1, c2);
  }

  // recorrer los numeros entre el primer y ultimo valor de la cadena
  for (j = end_num - 1; j > init_num; j--) {
    tmp  = values[c][j];
    //printf("\tValue = %llu\n", tmp);
    
    // recorrer los elem pares de bits (elem letras del genoma) del numero en cuestion
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
    search_gen();
  }
  
  // recorrer los n ultimos pares de bits del primer numero (inicio de la cadena)
  tmp  = values[c][init_num];
  //printf("\tValue = %llu\n", tmp);
  for (i = elem; i >= init_num_pos; i--, pos++) {
    search_gen();
  }

  //printf("\n");  
  return;
}

//====================================================================================
//====================================================================================
//====================================================================================

int encode_context(char* filename, char* directory) {

  printf("Init Genome Compresion\n");

  FILE *f1, *f2, *f3, *f4;
  unsigned long long value, value2;
  size_t size[50] = {0};
  size_t size2;
  int chromosome, i, cont, cont2;
  char *line1, *line2;
  char *tmp;
  size_t contador = 100000000;

  int elem = sizeof(unsigned long long) << 2;
  //printf("Elem = %i\n", elem);

  size2 = strlen(directory);
  tmp = (char *)malloc((size2 + 40) * sizeof(char));
  line1 = (char *)malloc(512 * sizeof(char));
  line2 = (char *)malloc(512 * sizeof(char));

  f1 = fopen (filename, "r");
  if (f1==NULL) {
    perror("No se puede abrir el fichero de entrada");
    return -1;
  }
  sprintf(tmp, "%s/Genome_context_CT.bin", directory);
  //printf("Fichero salida 1: %s\n", tmp);
  f2 = fopen (tmp, "wb");
  if (f2==NULL) {
    perror("No se puede abrir el fichero de contexto CT");
    return -1;
  }
  sprintf(tmp, "%s/Genome_context_GA.bin", directory);
  //printf("Fichero salida 2: %s\n", tmp);
  f3 = fopen (tmp, "wb");
  if (f3==NULL) {
    perror("No se puede abrir el fichero de contexto GA");
    return -1;
  }
  sprintf(tmp, "%s/Genome_context_size.txt", directory);
  //printf("Fichero salida 3: %s\n", tmp);
  f4 = fopen (tmp, "w");
  if (f4==NULL) {
    perror("No se puede abrir el fichero de tamaños");
    return -1;
  }
  
  chromosome = 0;
  size[chromosome] = 0;
  cont  = 0;
  cont2 = 0;

  // descartar la primera linea
  do {
    fgets(line1, 512, f1);
  } while(line1[0] == '>');
  
  
  while (fgets(line2, 512, f1) != NULL && contador) {
    contador--;
    if (line1[0] == '>') {
      size[chromosome]++;
      printf("chromosomes = %2i, size = %10lu\n", chromosome, size[chromosome]);
      chromosome++;
      for (; cont > 0 && cont <= elem; cont++) {
	  value = value << 2;
      }
      for (; cont2 > 0 && cont2 <= elem; cont2++) {
	  value2 = value2 << 2;
      }
      fwrite(&value,  sizeof(unsigned long long), 1, f2);
      fwrite(&value2, sizeof(unsigned long long), 1, f3);
      cont  = 0;
      cont2 = 0;
    } else {
      size2 = strlen(line1);
      //size[chromosome] += size2;

      // value for C->T conversion
      for (i = 0; i < size2 - 2; i++, cont++) {
	if (line1[i] == 'C') {
	  if (line1[i + 1] == 'G') {
	    value += 1;
	  } else {
	    if (line1[i + 2] == 'G') {
	      value += 2;
	    } else {
	      value += 3;
	    }
	  }
	}
	if (cont == elem) {
	  cont = 0;
	  size[chromosome]++;
	  fwrite(&value, sizeof(unsigned long long), 1, f2);
	} else {
	  value = value << 2;
	}
      }
      
      if (line1[size2 - 2] == 'C') {
	if (line1[size2 - 1] == 'G') {
	  value += 1;
	} else {
	  if (line2[0] == 'G') {
	    value += 2;
	  } else {
	    value += 3;
	  }
	}
      }
      if (cont == elem) {
	cont = 0;
	size[chromosome]++;
	fwrite(&value, sizeof(unsigned long long), 1, f2);
      } else {
	  value = value << 2;
      }
      cont++;

      if (line1[size2 - 1] == 'C') {
	if (line2[0] == 'G') {
	  value += 1;
	} else {
	  if (line2[1] == 'G') {
	    value += 2;
	  } else {
	    value += 3;
	  }
	}
      }
      if (cont == elem) {
	cont = 0;
	size[chromosome]++;
	fwrite(&value, sizeof(unsigned long long), 1, f2);
      } else {
	  value = value << 2;
      }
      cont++;
      // end value for C->T conversion

      // value for G->A conversion
      for (i = 2; i < size2; i++, cont2++) {
	if (line1[i] == 'G') {
	  if (line1[i - 1] == 'C') {
	    value2 += 1;
	  } else {
	    if (line1[i - 2] == 'C') {
	      value2 += 2;
	    } else {
	      value2 += 3;
	    }
	  }
	}
	if (cont2 == elem) {
	  cont2 = 0;
	  fwrite(&value2, sizeof(unsigned long long), 1, f3);
	} else {
	  value2 = value2 << 2;
	}
      }
      
      if (line2[0] == 'G') {
	if (line1[size2 - 1] == 'C') {
	  value2 += 1;
	} else {
	  if (line1[size2 - 2] == 'C') {
	    value2 += 2;
	  } else {
	    value2 += 3;
	  }
	}
      }
      if (cont2 == elem) {
	cont2 = 0;
	fwrite(&value2, sizeof(unsigned long long), 1, f3);
      } else {
	value2 = value2 << 2;
      }
      cont2++;

      if (line2[1] == 'G') {
	if (line2[0] == 'C') {
	  value2 += 1;
	} else {
	  if (line1[size2 - 1] == 'C') {
	    value2 += 2;
	  } else {
	    value2 += 3;
	  }
	}
      }
      if (cont2 == elem) {
	cont2 = 0;
	fwrite(&value2, sizeof(unsigned long long), 1, f3);
      } else {
	value2 = value2 << 2;
      }
      cont2++;
      // end value for G->A conversion

    }
    free(line1);
    line1 = strdup(line2);
  }

  //printf("chromosomes = %i\n", chromosome);  
  fprintf(f4, "%i\n", chromosome);
  for (i = 0; i < chromosome; i++) {
    fprintf(f4, "%lu\n", size[i]);
  }

  free(tmp);
  free(line1);
  free(line2);
  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f4);

  printf("End Genome Compresion\n");

  return 0;
}

//====================================================================================

int load_encode_context(char* directory, unsigned long long **valuesCT, unsigned long long **valuesGA) {

  //printf("Init Load Genome\n");

  FILE *f2, *f3, *f4;
  size_t size, size2 = strlen(directory);
  char *tmp = (char *)malloc((size2 + 40) * sizeof(char));
  int chromosome, i;

  sprintf(tmp, "%s/Genome_context_CT.bin", directory);
  //printf("Fichero salida 1: %s\n", tmp);
  f2 = fopen (tmp, "rb");
  if (f2==NULL) {
    perror("No se puede abrir el fichero de contexto CT");
    return -1;
  }
  sprintf(tmp, "%s/Genome_context_GA.bin", directory);
  f3 = fopen (tmp, "rb");
  if (f3==NULL) {
    perror("No se puede abrir el fichero de contexto GA");
    return -1;
  }
  sprintf(tmp, "%s/Genome_context_size.txt", directory);
  //printf("Fichero salida 3: %s\n", tmp);
  f4 = fopen (tmp, "r");
  if (f4==NULL) {
    perror("No se puede abrir el fichero de tamaños");
    return -1;
  }

  fscanf(f4, "%i\n", &chromosome);
  //printf("chromosomes %2i\n", chromosome);

  for (i = 0; i < chromosome; i++) {
    fscanf(f4, "%lu\n", &size);
    //printf("\tchromosome %2i (%10lu)\n", i, size);

    valuesCT[i] = (unsigned long long *)calloc(size, sizeof(unsigned long long));
    fread (valuesCT[i], sizeof(unsigned long long), size, f2);

    valuesGA[i] = (unsigned long long *)calloc(size, sizeof(unsigned long long));
    fread (valuesGA[i], sizeof(unsigned long long), size, f3);
  }

  free(tmp);
  fclose(f2);
  fclose(f3);
  fclose(f4);

  //printf("--------\nCT = %llu\n--------\nGA = %llu\n", valuesCT[0][100000], valuesGA[0][100000]);
  //printf("End Load Genome\n");

  return chromosome;
}

//====================================================================================


