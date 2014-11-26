/*
 * utils.c
 *
 *  Created on: Oct 30, 2014
 *      Author: sgallego
 */

#include "utils.h"

//----------------------------------------------------------------------
// reverse and convert secuence
//----------------------------------------------------------------------


void revcomp_seq(char* seq) {
  extern char convert_ASCII[128];

  char temp;
  int i, j = 0;
  
  i = 0;
  j = strlen(seq) - 1;
  
  while (i < j) {
    temp = seq[i];
    seq[i] = convert_ASCII[(int)seq[j]];
    seq[j] = convert_ASCII[(int)temp];
    i++;
    j--;
  }
}

//----------------------------------------------------------------------
// check if the file is pair
//----------------------------------------------------------------------
int is_pair(char *bam_filename){
	int p =0; //boolean que dice si es pair o no
	bam_file_t *bam_file = bam_fopen(bam_filename);
	bam1_t *bam1;
 	bam1 = bam_init1();
 	int i = 0;
 	while (!p && (bam_read1(bam_file->bam_fd, bam1) > 0) && (i<100)){
	//if (bam_read1(bam_file->bam_fd, bam1) > 0) {
	     uint32_t bam_flag = (uint32_t) bam1->core.flag;
	     i++;
	     if((!(bam_flag & BAM_FUNMAP)) && (bam_flag & BAM_FPAIRED)){
			 p = 1;
		 }

	}

	bam_destroy1(bam1);
	bam_fclose (bam_file);
	return p;
}

//------------------------------------------------------------------------
// check the file and return the max_quality
//------------------------------------------------------------------------

int max_quality(char *bam_filename){
	int q = 0; //máxima calidad
	bam_file_t *bam_file = bam_fopen(bam_filename);
	bam1_t *bam1;
 	bam1 = bam_init1();

	while (bam_read1(bam_file->bam_fd, bam1) > 0) {
	     uint32_t bam_flag = (uint32_t) bam1->core.flag;
	     if((!(bam_flag & BAM_FUNMAP))){
			 if(bam1->core.qual > q){
				 q = bam1->core.qual;
			 }
		 }
	}
	bam_destroy1(bam1);
	bam_fclose (bam_file);
	return q;
}

//------------------------------------------------------------------------

char *bam1_get_sequence(bam1_t *bam1, char *sequence_string) {
  /******************************************************
   * This function return the sequence in ASCII
   * @param bam1 bam for get the sequence
   * @param sequence_string secuencia previamente declarada; remember keep memory before call the function
   * *****************************************************/
  uint8_t* sequence_p = bam1_seq(bam1); //sam_tools: each base is encoded in 4 bits
  int sequence_length = (int32_t)bam1->core.l_qseq;
  unsigned char ta_lut[16];
  for (int i = 0; i < 16; i++){
    ta_lut [i] = 0;
  }
  ta_lut[1] = 'A';
  ta_lut[2] = 'C';
  ta_lut[4] = 'G';
  ta_lut[8] = 'T';
  ta_lut[15] = 'N';
  for (int i = 0; i < sequence_length; i++) {
    sequence_string[i] = ta_lut[bam1_seqi(sequence_p, i)];
  }
  sequence_string[sequence_length] = '\0';

  return sequence_string;
}

//------------------------------------------------------------------------

char *bam1_get_quality(bam1_t *bam1, char *quality_string){
  /********************************************************
   * This function convert the base Quality+33 in ASCII
   * @param bam1 bam for computer the sequence
   * @param quality_string remember keep memory before call the function
   * ******************************************************/
  int quality_length = (int32_t)bam1->core.l_qseq;
  for (int i = 0; i < quality_length; i++){
    quality_string [i] = bam1_qual(bam1)[i] + 33;
  }

  return quality_string;
}

//------------------------------------------------------------------------
