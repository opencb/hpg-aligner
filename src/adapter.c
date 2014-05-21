#include "adapter.h"
#include "bioformats/fastq/fastq_read.h"

#define ADAPTER_PREFIX 5

//--------------------------------------------------------------------

typedef struct adapter_match {
  int match;
  int adapter_start;
  int adapter_end;
  int seq_start;
  int seq_end;
  int num_mismatches;
  int length;
} adapter_match_t;

adapter_match_t *adapter_match_new() {
  adapter_match_t *p = (adapter_match_t *) calloc(1, sizeof(adapter_match_t));
  return p;
}

void adapter_match_free(adapter_match_t *p) {
  if (p) free(p);
}

void adapter_match_init(adapter_match_t *p) {
  //  assert(p != NULL);
  p->match = 0;
  p->num_mismatches = 0;
  p->length = 0;
}

void adapter_match_display(char *msg, adapter_match_t *p) {
  //  assert(p != NULL);
  if (p->match) {
    printf("%smatched (%i) at adapter [%i - %i] and sequence [%i - %i] (num. mismatches = %i)\n",
	   msg, p->length, p->adapter_start, p->adapter_end, p->seq_start, p->seq_end, p->num_mismatches);
  } else {
    printf("%sno matched (length = %i, num. mismatches = %i)\n", msg, p->length, p->num_mismatches);
  }
}

//--------------------------------------------------------------------

void match_adapter(char *adapter, int adapter_length, 
		   char *sequence, int sequence_length,
		   float ratio_error, adapter_match_t *match) {

  int i, j, len, pos, num_mismatches;
  char *p, *seq = sequence;
  char prefix[ADAPTER_PREFIX + 1];
  memcpy(prefix, adapter, ADAPTER_PREFIX);
  prefix[ADAPTER_PREFIX] = 0;
  //  printf("\t\tprefix = %s\n", prefix);

  match->match = 0;
  match->num_mismatches = 0;
  
  while ((p = strstr(seq, prefix)) != NULL) {
    pos = p - sequence;
    len = 0;
    num_mismatches = 0;
    for (i = pos + ADAPTER_PREFIX, j = ADAPTER_PREFIX; 
	 i < sequence_length && j < adapter_length; 
	 i++, j++) {
      if (sequence[i] != adapter[j]) {
	num_mismatches++;
      }
      len++;
    }
    match->length = ADAPTER_PREFIX + len;
    match->num_mismatches = num_mismatches;
    if (num_mismatches <= round(ratio_error * len)) {
      match->match = 1;
      match->adapter_start = 0;
      match->adapter_end = j - 1;
      match->seq_start = pos;
      match->seq_end = i - 1;
      break;
    }
    seq = p + 1;
  }
}

//--------------------------------------------------------------------

void cut_adapter(char *adapter, int adapter_length, fastq_read_t *read) {

  int len, pos;
  adapter_match_t *match = adapter_match_new();

  // search adapter in the 'forward' sequence
  //  printf("seq    : %s (adapter %s)\n", read->sequence, adapter);
  match_adapter(adapter, adapter_length, read->sequence, read->length, 0.1f, match);
  //  adapter_match_display("\t\t", match);
  if (match->match) {
    read->adapter_strand = 0;
    pos = read->length - (read->length / 3);
    fastq_read_display(read);
    //    printf("trimming adapter...\n");
    if (match->seq_start > pos || match->seq_end > pos) {

      // set adapter sequence and quality
      read->adapter_length = read->length - match->seq_start;
      read->adapter = (char *) malloc(read->adapter_length + 1);
      strcpy(read->adapter, &read->sequence[match->seq_start]);
      read->adapter[read->adapter_length] = 0;

      read->adapter_quality = (char *) malloc(read->adapter_length + 1);
      strcpy(read->adapter_quality, &read->quality[match->seq_start]);
      read->adapter_quality[read->adapter_length] = 0;

      // update sequence, quality and length
      read->sequence[match->seq_start] = 0;
      read->quality[match->seq_start] = 0;
      read->length -= read->adapter_length;

      // update reverse-complementary (adapter and sequence)
      read->adapter_revcomp = (char *) malloc(read->adapter_length + 1);
      strncpy(read->adapter_revcomp, read->revcomp, read->adapter_length);
      read->adapter_revcomp[read->adapter_length] = 0;

      strcpy(read->revcomp, &read->revcomp[read->adapter_length]);
      read->revcomp[read->length] = 0;

    } else {
      pos = read->length / 3;
      if (match->seq_start < pos || match->seq_end < pos) {
	// set adapter sequence and quality
	read->adapter_length = match->seq_end + 1;
	read->adapter = (char *) malloc(read->adapter_length + 1);
	strncpy(read->adapter, read->sequence, read->adapter_length);
	read->adapter[read->adapter_length] = 0;

	read->adapter_quality = (char *) malloc(read->adapter_length + 1);
	strncpy(read->adapter_quality, read->quality, read->adapter_length);
	read->adapter_quality[read->adapter_length] = 0;
	
	// update sequence, quality and length
	read->length -= read->adapter_length;
	strcpy(read->sequence, &read->sequence[read->adapter_length]);
	read->sequence[read->length] = 0;
	strcpy(read->quality, &read->quality[read->adapter_length]);
	read->quality[read->length] = 0;
	
	// update reverse-complementary (adapter and sequence)
	read->adapter_revcomp = (char *) malloc(read->adapter_length + 1);
	strcpy(read->adapter_revcomp, &read->revcomp[read->length]);
	read->adapter_revcomp[read->adapter_length] = 0;

	read->revcomp[read->length] = 0;

	read->adapter_length *= (-1);
      } else {
	//	printf("Not implemented yet!!\n");
	//	abort();
      }
    }
    //    fastq_read_display(read);      
  } else {
    // search adapter in the reverse-complementary sequence
    //    printf("revcomp: %s (adapter %s)\n", read->revcomp, adapter);
    match_adapter(adapter, adapter_length, read->revcomp, read->length, 0.1f, match);
    //    adapter_match_display("\t\t", match);
    if (match->match) {
      read->adapter_strand = 1;
      pos = read->length - (read->length / 3);
      //      fastq_read_display(read);
      //      printf("trimming adapter...\n");
      if (match->seq_start > pos || match->seq_end > pos) {
	
	// set adapter sequence and quality
	read->adapter_length = read->length - match->seq_start;
	read->adapter_revcomp = (char *) malloc(read->adapter_length + 1);
	strcpy(read->adapter_revcomp, &read->revcomp[match->seq_start]);
	read->adapter_revcomp[read->adapter_length] = 0;

	read->adapter_quality = (char *) malloc(read->adapter_length + 1);
	strcpy(read->adapter_quality, &read->quality[match->seq_start]);
	read->adapter_quality[read->adapter_length] = 0;
	
	// update sequence, quality and length
	read->revcomp[match->seq_start] = 0;
	read->length -= read->adapter_length;
	
	// update its reverse-complementary (adapter and sequence)
	read->adapter = (char *) malloc(read->adapter_length + 1);
	strncpy(read->adapter, read->sequence, read->adapter_length);
	read->adapter[read->adapter_length] = 0;

	strcpy(read->sequence, &read->sequence[read->adapter_length]);
	read->sequence[read->length] = 0;
	strcpy(read->quality, &read->quality[read->adapter_length]);
	read->quality[read->length] = 0;
	
      } else {
	pos = read->length / 3;
	if (match->seq_start < pos || match->seq_end < pos) {
	  // set adapter sequence and quality
	  read->adapter_length = match->seq_end + 1;
	  read->adapter_revcomp = (char *) malloc(read->adapter_length + 1);
	  strncpy(read->adapter_revcomp, read->revcomp, read->adapter_length);
	  read->adapter_revcomp[read->adapter_length] = 0;

	  read->adapter_quality = (char *) malloc(read->adapter_length + 1);
	  strncpy(read->adapter_quality, read->quality, read->adapter_length);
	  read->adapter_quality[read->adapter_length] = 0;
	  
	  // update sequence, quality and length
	  read->length -= read->adapter_length;
	  strcpy(read->revcomp, &read->revcomp[read->adapter_length]);
	  read->revcomp[read->length] = 0;
	  
	  // update its reverse-complementary (adapter and sequence)
	  read->adapter = (char *) malloc(read->adapter_length + 1);
	  strcpy(read->adapter, &read->sequence[read->length]);
	  read->adapter[read->adapter_length] = 0;

	  read->sequence[read->length] = 0;
	  read->quality[read->length] = 0;
	  
	  read->adapter_length *= (-1);
	} else {
	  //	  printf("Not implemented yet!!\n");
	  //	  abort();
	}
      }
      //      fastq_read_display(read);      
    } else {
      /*
      printf("*******************\n");
      printf("*******************\n");
      printf("*************** >>>>>>> no matched <<<<<<< *****************!!\n");
      printf("*******************\n");
      printf("*******************\n");
      */
    }
  }
  
  adapter_match_free(match);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
