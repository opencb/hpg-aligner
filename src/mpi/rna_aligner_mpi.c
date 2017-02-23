#include "rna_aligner_mpi.h"

#include <fcntl.h>

#define TMP_PATH_FILES    "TMP_FILES"
#define TMP_PATH_BUFFERS  "TMP_BUFFERS"

//#define TMP_PATH_OUTPUTS  "TMP_OUTPUTS"

#define MAX_BUFFER       10485760 //10MB of chars
//#define MAX_BUFFER       20971520 //20MB of chars
//#define MAX_BUFFER       31457280 //30MB of chars

//#define MAX_BUFFER         104857600 //100MB of chars

#define TOTAL_BUFFERS      3

#define MODE_1             0
#define MODE_2             1

#define SAME_NODE          0
#define SEND_TASK          1 
#define REQUEST_TASK       2
#define END_MPI            3
#define KILL_TH            4


#define MPI_TAG            1
#define MPI_TAG_MNG        2
#define MPI_TAG_SEEK       3

MPI_Request request_send[TOTAL_BUFFERS];
int id_send;

typedef struct n_tasks {
  int num_tasks;
  int first_id;
} n_tasks_t;
  
typedef struct work_package {
  int mode;
  int source;
  int target;
  char task[1024];
} work_package_t;

typedef struct MPI_data {
  volatile int block;
  volatile int run;
  pthread_mutex_t queue_mutex;
  pthread_mutex_t run_mutex;
  linked_list_t *queue_tasks;
} MPI_data_t;


typedef struct batch_out {
  char *data;
  size_t max_len;
  size_t len;
} batch_out_t;


typedef struct batch_buffer {
  array_list_t *buffer;
} batch_buffer_t;


int mpi_tot_send;
pthread_mutex_t mutex_merge;
int mpi_tot_reads;
double t_mpi_send;

size_t reads_wf;
size_t reads_r;

//====================================//
// R E A D E R    P A R A M E T E R S //
size_t start_seek;
size_t end_seek;
size_t r_tot_reads;
//====================================//

batch_out_t *new_batch_out() {
  batch_out_t *new_batch = (batch_out_t *)malloc(sizeof(batch_out_t));
  
  new_batch->data = (char *)calloc(MAX_BUFFER, sizeof(char));
  new_batch->len = 0;
  new_batch->max_len = MAX_BUFFER;
  
  return new_batch;
  
}

void free_batch_out(batch_out_t *batch_out) {
  free(batch_out->data);
  free(batch_out);
}

batch_buffer_t *new_batch_buffer() {
  batch_buffer_t *bb = (batch_buffer_t *)malloc(sizeof(batch_buffer_t));  
  bb->buffer = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  batch_out_t *b = new_batch_out();
  array_list_insert(b, bb->buffer);
  
  //for (int i = 0; i < TOTAL_BUFFERS; i++) {
  //b = new_batch_out();
  //array_list_insert(b, bb->buffer);
  //}
  
  b = new_batch_out();
  array_list_insert(b, bb->buffer);
  
  return bb;
  
}

void clear_batch_buffer(batch_buffer_t *bb) {
  int items = array_list_size(bb->buffer);
  for (int i = 0; i < items; i++) {
    batch_out_t *out = array_list_get(i, bb->buffer);
    free_batch_out(out);
  }
}

void free_batch_buffer(batch_buffer_t *bb) {
  array_list_free(bb->buffer, (void *)NULL);
  free(bb);
}


void write_to_buffer(batch_buffer_t *data_out, char *str) {  
  //printf("Consumer...\n");
  int len = strlen(str);
  int flag = 0;
  //int req_pos = id_send % TOTAL_BUFFERS;  
  //batch_out_t *b = array_list_get(req_pos + 1, data_out->buffer);
  batch_out_t *b = array_list_get(1, data_out->buffer);
  MPI_Status status;
  int rank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if ((b->len + len) > b->max_len) {
    int numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    //mpi_tot_send++;
    //struct timeval start, stop;
    //start_timer(start);

    //printf("[%i]Start Send id_send = %i\n", rank, id_send);

    MPI_Send(b->data, MAX_BUFFER, MPI_CHAR, numprocs - 1, 1, MPI_COMM_WORLD);

    //MPI_Isend(b->data, MAX_BUFFER, MPI_CHAR, numprocs - 1, 1, MPI_COMM_WORLD, &request_send[req_pos]);    
    //stop_timer(start, stop, t_mpi_send);    
    //MPI_Status status;
    //MPI_Wait(&request_send[req_pos], &status);    
    //printf("SEND and unlock!\n");
        
    b->len = 0;   
    /*
    id_send++;
    if (id_send >= TOTAL_BUFFERS) {
      req_pos = id_send % TOTAL_BUFFERS;
      MPI_Wait(&request_send[req_pos], &status);
      printf("[%i]Out wait.. req_pos = %i\n", rank, req_pos);
    }
    */
  }
  
  strncpy(&b->data[b->len], str, len);
  b->len += len;  
  
}


//--------------------------------------------------------------------


void bwt_convert_batch_to_str(pair_server_input_t *input, batch_t *batch, int max_alig) {

  //sw_server_input_t *sw_input = batch->sw_input;
  //genome_t *genome = sw_input->genome_p;
  //avls_list_t *avls_list = sw_input->avls_list;
  //metaexons_t *metaexons = sw_input->metaexons;
  
  batch_writer_input_t *writer_input = batch->writer_input;
  mapping_batch_t *mapping_batch = (mapping_batch_t *) batch->mapping_batch;
  int num_reads = array_list_size(batch->mapping_batch->fq_batch);
  //if (!num_reads) { return; }
  
  array_list_t *mapping_list;
  size_t num_mappings;
  
  pair_mng_t *pair_mng = input->pair_mng;
  report_optarg_t *report_optarg = input->report_optarg;
  int bam_format = batch->writer_input->bam_format;
  int pair_mode = pair_mng->pair_mode;  
  //size_t num_mappings;
  
  batch->data_output_size = 0;
  
  //num_reads = sa_batch->num_reads;
  
  int flag, pnext = 0, tlen = 0;
  char rnext[4] = "*\0";
  char optional_flags[512] = "\0";
  
  fastq_read_t *read;
  array_list_t *read_list = mapping_batch->fq_batch;
  alignment_t *alig;
  //array_list_t *mapping_list;

  sw_server_input_t *sw_input = batch->sw_input;
  genome_t *genome = sw_input->genome_p;
  
  //extern size_t total_reads, unmapped_reads, correct_reads;
  extern size_t total_reads;
  extern size_t reads_no_map;
  extern pthread_mutex_t mutex_sp;
  struct timeval time_free_s, time_free_e;
  extern double time_free_alig, time_free_list, time_free_batch;

  int num_mate_mappings;
  
  read = (fastq_read_t *) array_list_get(0, read_list);

  int read_size;
  int buffer_max_size = MAX_BUFFER;
  int buffer_size = 0;
  
  //Malloc Buffer
  //char *buffer = (char *)malloc(sizeof(char)*buffer_max_size);
  //buffer[0] = '\0';
    
  extern basic_statistics_t *basic_st;
  total_reads = num_reads;
  size_t num_mapped_reads = 0;
  size_t total_mappings = 0;
  size_t reads_uniq_mappings = 0;
  int buffer_pos;

  //Final
  batch_buffer_t *data_out = batch->data_out;
  //Temp
  //
  
  //mpi_tot_reads += num_reads;

  batch_out_t *b = array_list_get(0, data_out->buffer);
  
  char *buffer_tmp = b->data;
  size_t buffer_len = 0;

  reads_wf += num_reads;
  
  //printf("Start 1\n");
  for (size_t i = 0; i < num_reads; i++) {
    read = (fastq_read_t *) array_list_get(i, read_list);
    mapping_list = mapping_batch->mapping_lists[i];
    num_mappings = array_list_size(mapping_list);
    
    int mqual = 0;
    if (num_mappings == 1) { 
      reads_uniq_mappings++; 
      mqual = 50;
    } else if (num_mappings == 2) {
      mqual = 3;
    } else if (num_mappings == 3 || 
	       num_mappings == 4) {
      mqual = 1;
    }

    read_size   = ((strlen(read->id) + 50) + ((read->length + 50)*2) + (1000));
    
    if (num_mappings > 0) {
      num_mapped_reads++;
      total_mappings += num_mappings;
      for (size_t j = 0; j < num_mappings; j++) {
	alig = (alignment_t *) array_list_get(j, mapping_list);
	flag = (alig->seq_strand ? 16 : 0);

	sprintf(optional_flags, "AS:%i NM:%i NH:%lu", alig->map_quality, alig->mapq, num_mappings);

	read_size   = ((strlen(read->id) + 50) + ((read->length + 50)*2) + (1000));
	char alig_str[read_size];
	alig_str[0] = '\0';

	//printf("Alig %lu\n", alig->chromosome );
	sprintf(alig_str, "%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n%c", 
		alig->query_name,
		flag,
		genome->chr_name[alig->chromosome],
		alig->position + 1,
		mqual,
		alig->cigar,
		rnext,
		pnext,
		tlen,
		alig->seq_strand ? read->revcomp : read->sequence,
		//alig->sequence,
		alig->quality,
		optional_flags,
		'\0'
		);       
	
	strcpy(&buffer_tmp[buffer_len], alig_str);
	buffer_len += strlen(alig_str);
	assert(buffer_len < MAX_BUFFER);

	alignment_free(alig);
	
      }
    } else {

      char alig_str[read_size];

      alig_str[0] = '\0';      
      sprintf(alig_str, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n%c", 
	      read->id,
	      read->sequence,
	      read->quality,
	      '\0'
	      );

      strcpy(&buffer_tmp[buffer_len], alig_str);
      buffer_len += strlen(alig_str);
      assert(buffer_len < MAX_BUFFER);      

    }
    
    array_list_free(mapping_list, (void *) NULL);
    
  }
  //printf("1 - END\n");
  
  //batch->data_output = buffer;
  //batch->data_output_size = buffer_size;
  
  // free memory  
  if (mapping_batch) {
    mapping_batch_free(mapping_batch);
  }

  //Merge to Final buffer
  //pthread_mutex_lock(&mutex_merge);
  write_to_buffer(data_out, buffer_tmp);
  //pthread_mutex_unlock(&mutex_merge);

  //  free(buffer_tmp);
  
  //basic_statistics_add(total_reads, num_mapped_reads, total_mappings, reads_uniq_mappings, basic_st);

  return;

}


void MPI_output(void *data) {
  batch_t *batch = (batch_t *) data;

  //printf("CONVERT\n");
  bwt_convert_batch_to_str(batch->pair_input, batch, 100);  
  //free(batch->data_output);
  
  if (batch) batch_free(batch);

  //return CONSUMER_STAGE;
}

/*
void MPI_output(void *data) {
  batch_t *batch = (batch_t *) data;  
  batch_writer_input_t *writer_input = batch->writer_input;
  //bam_file_t *bam_file = writer_input->bam_file;    
  FILE *out_file = (FILE *) writer_input->bam_file;
  genome_t *genome = writer_input->genome;
  fastq_read_t *read;
  mapping_batch_t *mapping_batch = (mapping_batch_t *) batch->mapping_batch;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  //linked_list_t *linked_list = writer_input->list_p;
  array_list_t *read_list = mapping_batch->fq_batch;
  array_list_t *mapping_list;
  size_t num_mappings;
  size_t num_mapped_reads = 0;
  size_t total_mappings = 0;  
  alignment_t *alig;
  int flag, pnext = 0, tlen = 0;
  char rnext[4] = "*\0";
  
  
  for (size_t i = 0; i < num_reads; i++) {
    read = (fastq_read_t *) array_list_get(i, read_list);
    mapping_list = mapping_batch->mapping_lists[i];
    num_mappings = array_list_size(mapping_list);

    for (size_t j = 0; j < num_mappings; j++) {
      alig = (alignment_t *) array_list_get(j, mapping_list);
      alignment_free(alig);	 
    }

    if (mapping_list) {
      array_list_free(mapping_list, NULL);
    }	 

  }
  
  if (mapping_batch) {
    mapping_batch_free(mapping_batch);
  }

  if (batch) batch_free(batch);
  
  
  return CONSUMER_STAGE;

}
*/

size_t MPI_fastq_fread_bytes_se(array_list_t *reads, size_t bytes, fastq_file_t *fq_file) {
  size_t accumulated_size = 0;
  char header1[MAX_READ_ID_LENGTH];
  char sequence[MAX_READ_SEQUENCE_LENGTH];
  char header2[MAX_READ_ID_LENGTH];
  char qualities[MAX_READ_SEQUENCE_LENGTH];
  int header_length, sequence_length, quality_length, header2_length;
  fastq_read_t *read;
  char *res;

  
  while ((start_seek < end_seek) && (accumulated_size < bytes) &&
	 (fgets(header1, MAX_READ_ID_LENGTH, fq_file->fd) != NULL)) {
    
    res = fgets(sequence, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);
    res = fgets(header2, MAX_READ_ID_LENGTH, fq_file->fd);
    res = fgets(qualities, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);

    header_length   = strlen(header1);
    sequence_length = strlen(sequence);
    header2_length  = strlen(header2);
    quality_length  = strlen(qualities);

    // '\n' char is removed, but '\0' is left
    chomp_at(header1, header_length - 1);
    chomp_at(sequence, sequence_length - 1);
    chomp_at(qualities, quality_length - 1);
		
    read = fastq_read_new(header1, sequence, qualities);
    array_list_insert(read, reads);

    start_seek += header_length + sequence_length + header2_length + quality_length;
    
    accumulated_size += header_length + sequence_length + quality_length;
  }

  reads_r += array_list_size(reads);
  
  return accumulated_size;
  
}


void* SA_MPI_producer(void *input) {
  sa_wf_input_t *wf_input = (sa_wf_input_t *) input;
  sa_wf_batch_t *new_wf_batch = NULL;
  sa_wf_batch_t *curr_wf_batch = wf_input->wf_batch;  
  fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;
  array_list_t *reads = array_list_new(fq_reader_input->batch_size, 1.25f, 
				       COLLECTION_MODE_ASYNCHRONIZED);
  
  //extern size_t fd_read_bytes;
  
  // Fastq file
  if (fq_reader_input->flags == SINGLE_END_MODE) {
    MPI_fastq_fread_bytes_se(reads, fq_reader_input->batch_size, fq_reader_input->fq_file1);
  }  

  //printf("Total reads %lu\n", array_list_size(reads));
  
  size_t num_reads = array_list_size(reads);
  r_tot_reads += num_reads;
  
  if (num_reads == 0) {
    array_list_free(reads, (void *)fastq_read_free);
  } else {
    sa_batch_t *sa_batch = sa_batch_new(reads);
    
    new_wf_batch = sa_wf_batch_new(NULL,
				   curr_wf_batch->sa_index,
				   curr_wf_batch->writer_input, 
				   sa_batch,
				   curr_wf_batch->data_input);
    new_wf_batch->buffer_mpi = curr_wf_batch->buffer_mpi;
  }

  return new_wf_batch;

}


void* BWT_MPI_producer(void *input) {
  wf_input_t *wf_input = (wf_input_t *) input;
  batch_t *new_batch = NULL;
  batch_t *batch = wf_input->batch;
  fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;
  array_list_t *reads = array_list_new(10000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  
  MPI_fastq_fread_bytes_se(reads, fq_reader_input->batch_size, fq_reader_input->fq_file1);
  
  size_t num_reads = array_list_size(reads);
  
  if (num_reads == 0) {
    array_list_free(reads, (void *)fastq_read_free);
  } else {
    mapping_batch_t *mapping_batch = mapping_batch_new(reads, 
						       batch->pair_input->pair_mng);
    
    new_batch = batch_new(batch->bwt_input, batch->region_input, batch->cal_input, 
			  batch->pair_input, batch->preprocess_rna, batch->sw_input, batch->writer_input, 
			  batch->mapping_mode, mapping_batch, batch->data_out);
  }
  
  return new_batch;
  
}



void SA_MPI_consumer(void *data) {
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  char *str_output = (char *)wf_batch->data_output;
  size_t str_len = wf_batch->data_output_size;
  
  batch_buffer_t *buffer_mpi = wf_batch->buffer_mpi;
    
  if (!wf_batch->data_output_size) { return 0; }

  write_to_buffer(buffer_mpi, str_output);
  free(str_output);
  
  if (wf_batch) sa_wf_batch_free(wf_batch);

  return 0; 
}

void MPI_consumer(void *data) {

  return;

  batch_t *batch = (batch_t *)data;
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  pair_server_input_t *pair_input = batch->pair_input;
  //pair_mng_t *pair_mng = pair_input->pair_mng;
  //report_optarg_t *report_optarg = pair_input->report_optarg;
  int bam_format = batch->writer_input->bam_format;
  int pair_mode = pair_input->pair_mng->pair_mode;
  batch_buffer_t *data_out = batch->data_out;
  
  //size_t num_reads, num_mappings;
  //batch_buffer_t *data_out = batch->data_out;
  
  int flag, pnext = 0, tlen = 0;
  char rnext[4] = "*\0";
  char optional_flags[512] = "\0";
  
  fastq_read_t *read;
  array_list_t *read_list = mapping_batch->fq_batch;
  alignment_t *alig;
  array_list_t *mapping_list;
  
  genome_t *genome = batch->writer_input->genome;
  
  if (bam_format || (pair_mode != SINGLE_END_MODE)) {
    LOG_FATAL("NOT IMPLEMENTED YET!\n");
  }


  write_to_buffer(data_out, batch->data_output);
  free(batch->data_output);

  
  // free memory
  //if (mapping_batch)
  //mapping_batch_free(mapping_batch);
  //}

  if (batch) batch_free(batch);

  return 0;

}


array_list_t* new_steal_targets(int numprocs, int rank, int mode) {
  array_list_t *targets = array_list_new(numprocs, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  if (mode == MODE_1) {
    /**************************************************************
     *           All Ranks Start Stealing to First Rank (0)       *
     *                 [ 0 | 1 | 2 | 3 | 4 | 5 ]                  *
     *                   <---|---|---|---|---|                    *
     **************************************************************/
    for (int i = 0; i < numprocs - 1; i++) {
      if (i != rank) {
	printf("[%i]r-target %i (numprocs %i)\n", rank, i, numprocs);
	array_list_insert(i, targets);
      }
    }
  } else if (mode == MODE_2) {
    /**************************************************************
     *             All Ranks Start Stealing to next Rank (n+1)    *
     *                 [ 0 | 1 | 2 | 3 | 4 | 5 ]                  *
     *                 ->|-->|-->|-->|-->|-->|-                   *
     **************************************************************/
    for (int i = 0; i < numprocs - 1; i++) {
      int target = (i + 1) % (numprocs - 1);
      if (target != rank) {
	printf("[%i]r-target %i (numprocs %i)\n", rank, target, numprocs);
	array_list_insert(target, targets);
      }
    }
  }
  
  return targets;
  
}


void mpi_thread_polling(void *input) {
  //===============================//
  //= P O L L I N G - T H R E A D =//
  //===============================//
  MPI_data_t *MPI_data = (MPI_data_t *)input;
  
  int tag = 1;
  MPI_Status status;
  
  int numprocs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  array_list_t *targets = new_steal_targets(numprocs, rank, MODE_2);  
  int total_targets = array_list_size(targets);
  
  work_package_t work_request;
  int target_inx;
  int bad_request = 0;
  char *file;
  int delay = 2;

  //printf("START: MPI-Polling-Thread Rank-%i\n", rank);
  
  while (1) {
    MPI_Recv(&work_request, sizeof(work_package_t), MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
    //printf("[%i]Polling RECV from %i (%i)\n", rank, work_request.source, work_request.mode);
    //sleep(delay);
    if (work_request.mode == SAME_NODE) {
      //Thread Worker without work... Polling thread steal tasks
      //sleep(1);
      if (!array_list_size(targets)) {
	MPI_data->run = 0;
	MPI_data->block = 0;
	continue;
	//break;
      }
      work_request.mode = REQUEST_TASK;
      target_inx = 0;
      work_request.source = rank;
      work_request.target = array_list_get(target_inx, targets);
      MPI_Send(&work_request, sizeof(work_package_t), MPI_CHAR, work_request.target, 1, MPI_COMM_WORLD);
      printf("\t\tPolling.Rank %i -> Send to %i File request\n", rank, work_request.target);
      //sleep(delay);
    } else if (work_request.mode == REQUEST_TASK) {
      //Steal Request
      //sleep(1);
      pthread_mutex_lock(&(MPI_data->queue_mutex));
      if (MPI_data->queue_tasks->size > 1) {
	file = (char *)linked_list_remove_last(MPI_data->queue_tasks);
	strcpy(work_request.task, file);
      } else {
	printf("[%i] No tasks, Send\n", rank);
	work_request.task[0] = '\0';
      }
      pthread_mutex_unlock(&(MPI_data->queue_mutex));

      work_request.mode   = SEND_TASK;
      work_request.target = work_request.source;
      work_request.source = rank;
	  
      MPI_Send(&work_request, sizeof(work_package_t), MPI_CHAR, work_request.target, 1, MPI_COMM_WORLD);
      printf("\t\tPolling.Rank %i -> Send to %i file: %s\n", rank, work_request.target, work_request.task);
      //sleep(delay);
    } else if (work_request.mode == SEND_TASK) {
      //sleep(1);
      if (work_request.task[0] == '\0') {
	//Request to Other Target
	work_request.mode = REQUEST_TASK;
	target_inx++;
	printf("[%i] - target inx %i >= total targets %i\n", rank, target_inx, total_targets);
	if (target_inx >= total_targets) {
	  //End Work
	  printf("[%i]Finalize with total targets %i\n", rank, target_inx);
	  work_request.mode   = END_MPI;
	  //printf("\t\tPooling.Rank %i target_inx = %i >= total_targets = %i, Loop 0...%i\n", rank, target_inx, total_targets, numprocs - 1);
	  for (int i = 0; i < numprocs - 1; i++) {
	    if (i != rank) {
	      work_request.target = i;
	      work_request.source = rank;
	      //printf("\t\t-->Pooling.Rank %i Send to %i\n", rank, i);
	      MPI_Send(&work_request, sizeof(work_package_t), MPI_CHAR, work_request.target, 1, MPI_COMM_WORLD);
	    }
	  }
	  MPI_data->run = 0;
	  MPI_data->block = 0;
	  //continue;
	} else {
	  work_request.target = array_list_get(target_inx, targets);
	  printf("[%i]Send to %i New file request\n", rank, work_request.target);
	  MPI_Send(&work_request, sizeof(work_package_t), MPI_CHAR, work_request.target, 1, MPI_COMM_WORLD);
	  //sleep(delay);
	}
      } else {
	printf("[%i]Steal Job : %s\n", rank, work_request.task);
	pthread_mutex_lock(&(MPI_data->queue_mutex));
	linked_list_insert_first(strdup(work_request.task), MPI_data->queue_tasks);
	pthread_mutex_unlock(&(MPI_data->queue_mutex));
	MPI_data->block = 0;
      }
    } else if (work_request.mode == END_MPI) {
      //sleep(1);
      printf("\t\tPolling.Rank %i END Loop\n", rank);
      MPI_data->run = 0;
      MPI_data->block = 0;
      continue;
      //break;
    } else {
      break;
    }
  }

  //MPI_data->run = 0;
  //MPI_data->block = 0;
  
  printf("\t[%i]END MPI-Polling-Thread\n", rank);

  
}



//Return TMP_PATH name
int split_input_file(char *fq_str, char *tmp_path, int reads_split) {

  FILE *tmp_fd;
  FILE *fq_file = fopen(fq_str, "r");

  //Prepare PATH tmp files
  create_directory(tmp_path);
  
  int count;
  char header1[MAX_READ_ID_LENGTH];
  char sequence[MAX_READ_SEQUENCE_LENGTH];
  char header2[MAX_READ_ID_LENGTH];
  char qualities[MAX_READ_SEQUENCE_LENGTH];

  char *res;

  char tmp_file[1024];
  int file_id = 0;
  res = fgets(header1, MAX_READ_ID_LENGTH, fq_file);
  
  while (res) {    
    sprintf(tmp_file, "%s/%i.tmp", tmp_path, file_id);
    tmp_fd = fopen(tmp_file, "w");
    count = 0;
    
    while (count < reads_split) {      
      res = fgets(sequence, MAX_READ_SEQUENCE_LENGTH, fq_file);
      res = fgets(header2, MAX_READ_ID_LENGTH, fq_file);
      res = fgets(qualities, MAX_READ_SEQUENCE_LENGTH, fq_file);

      fputs(header1, tmp_fd);
      fputs(sequence, tmp_fd);
      fputs(header2, tmp_fd);
      fputs(qualities, tmp_fd);
      
      count++;

      res = fgets(header1, MAX_READ_ID_LENGTH, fq_file);
      if (!res) break;
    }
    
    fclose(tmp_fd);   
    //if (!res) break;    
    file_id++;    
    
  }
  
  
  
  fclose(fq_file);

  //file_id++;
  //printf("Split file in %i\n", file_id);
  
  return file_id;
  
}

void *write_sam_header_BWT(genome_t *genome, FILE *f);


    /*
    unsigned long n_items = 1000000;
    unsigned long *test1 = (unsigned long *)malloc(sizeof(unsigned long)*n_items);
    unsigned long *test2 = (unsigned long *)malloc(sizeof(unsigned long)*n_items*10);
	
    unsigned long* recv_items, *recv_items2, *recv_items3;
	
    for (int t = 0; t < n_items; t++) {
      test1[t] = t;
    }
	
    for (int t = 0; t < n_items*10; t++) {
      test2[t] = t;
    }
	
    if (rank == 0) {
      MPI_Send(test1, sizeof(unsigned long)*n_items, MPI_BYTE, 1, tag, MPI_COMM_WORLD);
      MPI_Send(test2, sizeof(unsigned long)*n_items*10, MPI_BYTE, 1, tag, MPI_COMM_WORLD);
    } else {
      recv_items = (unsigned long *)malloc(sizeof(unsigned long)*n_items);
      MPI_Recv(recv_items, sizeof(unsigned long)*n_items, MPI_BYTE, 0, tag, MPI_COMM_WORLD, &status); 
      for (int t = 0; t < n_items; t++) {
	if (recv_items[t] != test1[t]) { printf("ERROR in Send 1\n"); exit(-1); }
      }
      printf("TEST1: OK\n");
	
      recv_items2 = (unsigned long *)malloc(sizeof(unsigned long)*n_items*10);
      MPI_Recv(recv_items2, sizeof(unsigned long)*n_items*10, MPI_BYTE, 0, tag, MPI_COMM_WORLD, &status); 
      for (int t = 0; t < n_items*10; t++) {
	if (recv_items2[t] != test2[t]) { printf("ERROR in Send 2\n"); exit(-1); }
      }
      printf("TEST2: OK\n");
	
      for (int t = 0; t < n_items*10; t++) {
	test2[t] = 0;
      }
	
    }
	
    MPI_Bcast(test2, sizeof(unsigned long)*n_items*10, MPI_BYTE, 0, new_comm);    
	
    if (rank != 0) {
      for (int t = 0; t < n_items*10; t++) {
	if (recv_items2[t] != test2[t]) { printf("ERROR in Send 3\n"); exit(-1); }
      }
      printf("TEST3: OK\n");      
    }    
    */
    //exit(-1);

/*
    //==== W O R K E R    E N D    P R O C E S S    F I L E ====//
    //====             S T A R T   M E R G E                ====//
    //==========================================================//
    
    //==========================================================================//
    //= M E R G E    M E T A E X O N - A V L     T O    2nd    W O R K F L O W =//
    //==========================================================================//
    //printf("[%i]Start merge...\n", rank);
    start_timer(wx_start);

    unsigned long num_sj;        
    int numprocs_merge = numprocs - 1;
    int dec_s = 1;
    int dec_r = 2;
    int send, recv;
    int dec = 2;
      
    while (1) {
      send = numprocs_merge - dec_s;
      recv = numprocs_merge - dec_r;
      if (send > 0 && recv < 0)
	recv = 0;
      if (send <= 0 && recv <= 0)
	break;
      
      for (int r = numprocs_merge - 1; r >= 0; r--) {
	if (rank == send) {
	  //SEND AVL
	  num_sj = get_total_items(num_chromosomes, avls_list);
	  unsigned long  *list_sj = (unsigned long *)malloc(sizeof(unsigned long)*num_sj*5); 	  
	  MPI_avl_package(num_chromosomes, avls_list, num_sj, list_sj);	  
	  MPI_Send(&num_sj, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);	  
	  
	  MPI_Send(list_sj, num_sj*5, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  
	  //SEND METAEXON
	  unsigned long num_meta = 0;
	  unsigned long num_left_breaks = 0, num_right_breaks = 0;
	  for (int c = 0; c < metaexons->num_chromosomes; c++) {
	    linked_list_t *metaexon_list = metaexons->metaexons_list[c];
	    linked_list_item_t *item;
	    for (item = metaexon_list->first; item != NULL; item = item->next) {
	      metaexon_t *metaexon = item->item;
	      num_left_breaks  += array_list_size(metaexon->left_breaks);
	      num_right_breaks += array_list_size(metaexon->right_breaks);
	      num_meta++;	      
	    }
	  }

	  printf("Pack meta %i, %i, %i...\n", num_meta, num_left_breaks, num_right_breaks);

	  unsigned long *list_meta = (unsigned long *)malloc(sizeof(unsigned long)*num_meta*5);
	  unsigned long *list_left_breaks  = (unsigned long *)malloc(sizeof(unsigned long)*num_left_breaks*2);
	  unsigned long *list_right_breaks = (unsigned long *)malloc(sizeof(unsigned long)*num_right_breaks*2);

	  MPI_metaexon_package(metaexons, num_meta, list_meta, 
			       num_left_breaks, list_left_breaks, 
			       num_right_breaks, list_right_breaks);

	  printf("RECV-RANK0: %i, %i, %i, %i, %i\n", list_meta[num_meta], list_meta[num_meta + 1], list_meta[num_meta + 2], list_meta[num_meta + 3], list_meta[num_meta + 4]);
	  MPI_Send(&num_meta, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(list_meta, num_meta*5, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  
	  MPI_Send(&num_left_breaks, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(list_left_breaks, num_left_breaks*2, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  
	  //printf("SEND: %i, %i, %i, %i, %i\n", list_right_breaks[0], list_right_breaks[1], list_right_breaks[2], list_right_breaks[3], list_right_breaks[4]);
	  MPI_Send(&num_right_breaks, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(list_right_breaks, num_right_breaks*2, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);

	  free(list_sj);
	  free(list_meta);
	  free(list_left_breaks);
	  free(list_right_breaks);

	}
	
	if (rank == recv) {
	  //Recv AVLs
	  unsigned long recv_num_sj;	    
	  MPI_Recv(&recv_num_sj, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  printf("recv sj = %lu\n", recv_num_sj);

	  unsigned long *recv_list_sj = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_sj*5);
	  MPI_Recv(recv_list_sj, recv_num_sj*5, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);	  
	  
	  unsigned long recv_num_meta;
	  MPI_Recv(&recv_num_meta, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  unsigned long *recv_list_meta = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_meta*5);
	  MPI_Recv(recv_list_meta, recv_num_meta*5, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  
	  unsigned long recv_num_left_breaks;
	  MPI_Recv(&recv_num_left_breaks, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  unsigned long *recv_list_left_breaks = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_left_breaks*2);
	  MPI_Recv(recv_list_left_breaks, recv_num_left_breaks*2, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  
	  unsigned long recv_num_right_breaks;
	  MPI_Recv(&recv_num_right_breaks, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  unsigned long *recv_list_right_breaks = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_right_breaks*2);
	  MPI_Recv(recv_list_right_breaks, recv_num_right_breaks*2, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);

	  printf("RECV: %i, %i, %i, %i, %i\n", recv_list_meta[recv_num_meta], recv_list_meta[recv_num_meta + 1], recv_list_meta[recv_num_meta + 2], recv_list_meta[recv_num_meta + 3], recv_list_meta[recv_num_meta + 4]);	  


	  //printf("[%i]Merge %i sj:\n", rank, recv_num_sj);	  

	  //Merge AVLs
	  avl_node_t *node_avl_start, *node_avl_end;
	  MPI_splice_t MPI_sj;
	  for (size_t sj = 0; sj < recv_num_sj; sj++) {
	    MPI_sj.start          = recv_list_sj[sj];
	    MPI_sj.end            = recv_list_sj[recv_num_sj + sj];
	    MPI_sj.reads_number   = recv_list_sj[recv_num_sj*2 + sj];
	    MPI_sj.strand         = recv_list_sj[recv_num_sj*3 + sj];
	    MPI_sj.chr            = recv_list_sj[recv_num_sj*4 + sj];

	    //printf("\t[%i] %i/%i Merge (%i) %i: %lu - %lu\n", rank, sj, recv_num_sj, MPI_sj.strand, MPI_sj.chr, MPI_sj.start, MPI_sj.end);
	    
	    allocate_start_node(MPI_sj.chr - 1,
				MPI_sj.strand,
				MPI_sj.start,
				MPI_sj.end,
				MPI_sj.start,
				MPI_sj.end,
				FROM_READ,
				1,
				"",
				&node_avl_start,
				&node_avl_end,
				avls_list);
	    //Refresh MPI_sj.reads_number in nodes
	    //end_data_t *end_data = node_avl_start->data;	    
	  }
	  //Merge Metaexon
	  //sleep(1);
	  
	  //printf("[%i]Merge %i META:\n", rank, recv_num_meta);
	  size_t l_pos = 0, r_pos = 0;
	  MPI_metaexon_t mpi_meta;
	  for (size_t m = 0; m < recv_num_meta; m++) {
	    mpi_meta.chromosome = recv_list_meta[m];
	    mpi_meta.start      = recv_list_meta[recv_num_meta   + m];
	    mpi_meta.end        = recv_list_meta[recv_num_meta*2 + m];
	    mpi_meta.n_starts   = recv_list_meta[recv_num_meta*3 + m];
	    mpi_meta.n_ends     = recv_list_meta[recv_num_meta*4 + m];

	    //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
	    MPI_breaks_t bk;
	    for (size_t s = 0; s < mpi_meta.n_starts; s++) {
	      bk.pos    = recv_list_right_breaks[r_pos];
	      bk.strand = recv_list_right_breaks[recv_num_right_breaks + r_pos];
	      r_pos++;
	      //MPI_breaks_t bk = MPI_breaks[se_pos++];
	      //printf("\tM-STARTS-(%i)%lu\n", bk.strand, bk.pos);
	      cp_avltree *avl = avls_list->avls[bk.strand][mpi_meta.chromosome].avl;
	      avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	      assert(avl_node);
	      //if (!avl_node) {
	      //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
	      //printf("%i:%i:%i\n", m, r_pos - 1, bk.pos);
	      //exit(-1);
	      //}
	      metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			      METAEXON_RIGHT_END, avl_node, metaexons);
	    }
	    for (size_t s = 0; s < mpi_meta.n_ends; s++) {
	      bk.pos    = recv_list_left_breaks[l_pos];
	      bk.strand = recv_list_left_breaks[recv_num_left_breaks + l_pos];
	      l_pos++;
	      //MPI_breaks_t bk = MPI_breaks[se_pos++];
	      //printf("\t[%i]M-ENDS-(%i)%lu\n", rank, bk.strand, bk.pos);
	      cp_avltree *avl = avls_list->ends_avls[bk.strand][mpi_meta.chromosome].avl;
	      avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	      assert(avl_node);
	      metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			      METAEXON_LEFT_END, avl_node, metaexons);
	    }
	  }

	  free(recv_list_sj);
	  free(recv_list_meta);
	  free(recv_list_left_breaks);
	  free(recv_list_right_breaks);
 	  
	}
	send -= dec;
	recv -= dec;
	if (send <= 0) { send = -1; }
      }
      dec   *= 2;
      dec_s *= 2;
      dec_r *= 2;
    }
  
    
    printf("[%i]Merge 1 end\n", rank);
    
    //exit(-1);
    
    //int num_sj;
    unsigned long recv_num_sj, recv_n_metaexons, recv_n_starts_ends;
    //MPI_splice_t *MPI_splice;
    //MPI_metaexon_package_t *mpi_meta_pack;
    //MPI_metaexon_t *MPI_metaexon;
    //MPI_breaks_t *MPI_breaks;
    
    unsigned long *list_sj;
    unsigned long *list_meta;
    unsigned long *list_left_breaks;
    unsigned long *list_right_breaks;
    
    unsigned long num_meta = 0;
    unsigned long num_left_breaks = 0, num_right_breaks = 0;    

    if (rank == 0) {
      //PACK AVL
      num_sj = get_total_items(num_chromosomes, avls_list);
      list_sj = (unsigned long *)malloc(sizeof(unsigned long)*num_sj*5); 	  
      MPI_avl_package(num_chromosomes, avls_list, num_sj, list_sj);	  
      
      //PACK METAEXON
      for (int c = 0; c < metaexons->num_chromosomes; c++) {
	linked_list_t *metaexon_list = metaexons->metaexons_list[c];
	linked_list_item_t *item;
	for (item = metaexon_list->first; item != NULL; item = item->next) {
	  metaexon_t *metaexon = item->item;
	  num_left_breaks  += array_list_size(metaexon->left_breaks);
	  num_right_breaks += array_list_size(metaexon->right_breaks);
	  num_meta++;	      
	}
      }
      
      list_meta = (unsigned long *)malloc(sizeof(unsigned long)*num_meta*5);
      list_left_breaks  = (unsigned long *)malloc(sizeof(unsigned long)*num_left_breaks*2);
      list_right_breaks = (unsigned long *)malloc(sizeof(unsigned long)*num_right_breaks*2);

      MPI_metaexon_package(metaexons, num_meta, list_meta, 
			   num_left_breaks, list_left_breaks, 
			   num_right_breaks, list_right_breaks);


      //printf("RECV2-RANK0: %i, %i, %i, %i, %i\n", list_meta[num_meta], list_meta[num_meta + 1], list_meta[num_meta + 2], list_meta[num_meta + 3], list_meta[num_meta + 4]);

    }
    
    //printf("[%i]Bcast Send\n", rank);
    
    MPI_Bcast(&num_sj, 1, MPI_UNSIGNED_LONG, 0, new_comm);
    MPI_Bcast(&num_meta, 1, MPI_UNSIGNED_LONG, 0, new_comm);
    MPI_Bcast(&num_left_breaks, 1, MPI_UNSIGNED_LONG, 0, new_comm);
    MPI_Bcast(&num_right_breaks, 1, MPI_UNSIGNED_LONG, 0, new_comm);
    printf("[%i]Pack meta %i %i, %i, %i...\n", rank, num_sj, num_meta, num_left_breaks, num_right_breaks);
    
    //size_t xn = 100000000;
    size_t xn = num_sj*5 + num_meta*5 + num_left_breaks*2 + num_right_breaks*2;
    unsigned long *pack_send = (unsigned long *)malloc(sizeof(unsigned long)*xn);

    if (rank == 0) {
      memcpy(pack_send, list_sj, sizeof(unsigned long)*num_sj*5);
      memcpy(&pack_send[num_sj*5], list_meta, sizeof(unsigned long)*num_meta*5);
      memcpy(&pack_send[num_sj*5 + num_meta*5], list_left_breaks, sizeof(unsigned long)*num_left_breaks*2);
      memcpy(&pack_send[num_sj*5 + num_meta*5 + num_left_breaks*2], list_right_breaks, sizeof(unsigned long)*num_right_breaks*2);

      free(list_sj);
      free(list_meta);
      free(list_left_breaks);
      free(list_right_breaks);
    }

    MPI_Bcast(pack_send, xn, MPI_UNSIGNED_LONG, 0, new_comm);

    printf("[%i] RECV2: %lu, %lu, %lu \n", rank, pack_send[0], pack_send[1], pack_send[num_meta]);

    list_sj           = pack_send;
    list_meta         = &pack_send[num_sj*5];
    list_left_breaks  = &pack_send[num_sj*5 + num_meta*5];
    list_right_breaks = &pack_send[num_sj*5 + num_meta*5 + num_left_breaks*2];
  
      //  exit(-1);

    printf("[%i]Merge 2 %i:\n", rank, num_sj);

    //Merge AVLs
    avl_node_t *node_avl_start, *node_avl_end;
    MPI_splice_t MPI_sj;

    for (size_t sj = 0; sj < num_sj; sj++) {
      MPI_sj.start          = list_sj[sj];
      MPI_sj.end            = list_sj[num_sj + sj];
      MPI_sj.reads_number   = list_sj[num_sj*2 + sj];
      MPI_sj.strand         = list_sj[num_sj*3 + sj];
      MPI_sj.chr            = list_sj[num_sj*4 + sj];

      //printf("[%i]Merge 2 (%i) %i: %lu - %lu\n", rank, MPI_sj.strand, MPI_sj.chr, MPI_sj.start, MPI_sj.end);
      allocate_start_node(MPI_sj.chr - 1,
			  MPI_sj.strand,
			  MPI_sj.start,
			  MPI_sj.end,
			  MPI_sj.start,
			  MPI_sj.end,
			  FROM_READ,
			  1,
			  "",
			  &node_avl_start,
			  &node_avl_end,
			  avls_list);
      //Refresh MPI_sj.reads_number in nodes
      end_data_t *end_data = node_avl_start->data;	    
    }
    
    printf("[%i]Merge Meta %i\n", rank, num_meta);
    
    //Merge Metaexon
    size_t se_pos = 0;
    size_t r_pos = 0, l_pos = 0;
    MPI_metaexon_t mpi_meta;
    
    for (size_t m = 0; m < num_meta; m++) {
      mpi_meta.chromosome = list_meta[m];
      mpi_meta.start      = list_meta[num_meta   + m];
      mpi_meta.end        = list_meta[num_meta*2 + m];
      mpi_meta.n_starts   = list_meta[num_meta*3 + m];
      mpi_meta.n_ends     = list_meta[num_meta*4 + m];
      
      //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
      MPI_breaks_t bk;
      for (int s = 0; s < mpi_meta.n_starts; s++) {
	bk.pos    = list_right_breaks[r_pos];
	bk.strand = list_right_breaks[num_right_breaks + r_pos];
	r_pos++;

	//printf("\tM-STARTS-(%i)%lu\n", bk.strand, bk.pos);
	cp_avltree *avl = avls_list->avls[bk.strand][mpi_meta.chromosome].avl;
	avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	//assert(avl_node);
	if (!avl_node) {
	  printf("ERROR: (%i|%i)%lu\n", m, s, bk.pos);
	  exit(-1);
	}
	metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			METAEXON_RIGHT_END, avl_node, metaexons);
      }
      for (int s = 0; s < mpi_meta.n_ends; s++) {
	bk.pos    = list_left_breaks[l_pos];
	bk.strand = list_left_breaks[num_left_breaks + l_pos];
	l_pos++;

	//printf("\t[%i]M-ENDS-(%i)%lu\n", rank, bk.strand, bk.pos);
	cp_avltree *avl = avls_list->ends_avls[bk.strand][mpi_meta.chromosome].avl;
	avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	assert(avl_node);
	metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			METAEXON_LEFT_END, avl_node, metaexons);
      }
    }
    
    printf("Frees merge\n");
    free(pack_send);

  
    //====                E N D    M E R G E                ====//
    //==========================================================//
*/



void manager_function(char *file_name) {
  
  //MALLEABLE VARIABLE FOR SPLIT FILE
  int tag_mng  = MPI_TAG_MNG;
  int tag_seek = MPI_TAG_SEEK;

  size_t split_factor = 10; // Configurable (input parameter?)
  size_t read_positions[split_factor + 1];
  
  //===================================================================//
  //                               O M P                               //
  //===================================================================//
  //                S P L I T    I N P U T    F I L E                  //
  //===================================================================//
  
  FILE *fd = fopen(file_name, "r");
  if (!fd) { printf("ERROR opening file\n"); exit(-1); }
    
  fseek(fd, 0L, SEEK_END);
  size_t fd_total_bytes = ftell(fd);
  fseek(fd, 0L, SEEK_SET);
    
  int i;
  size_t seek_bytes = 0;
  
  size_t inc_bytes = fd_total_bytes / split_factor;
    
  size_t offset;
  char line[2048];
    
  for (i = 0; i < split_factor; i++) {
    offset = 0;
    fseek(fd, seek_bytes, SEEK_SET);
    fgets(line, 2048, fd);
    while (line[0] != '@') {
      offset += strlen(line);
      fgets(line, 2048, fd);
    }      
    seek_bytes += offset;
    read_positions[i] = seek_bytes;
    seek_bytes += inc_bytes;
  }
    
  read_positions[i] = fd_total_bytes;
  fclose(fd);
  
  //===================================================================//
  //                                                                   //
  //===================================================================//    

  size_t actual_tsk = 0;
  int petition;
  size_t s_positions[2];

  MPI_Status status;

  //printf("");
  while (actual_tsk < split_factor) {
    
    MPI_Recv(&petition, 1, MPI_INT, MPI_ANY_SOURCE, tag_mng, MPI_COMM_WORLD, &status);
    
    s_positions[0] = read_positions[actual_tsk];
    s_positions[1] = read_positions[actual_tsk + 1];
    
    MPI_Send(s_positions, 2*sizeof(size_t), MPI_CHAR, status.MPI_SOURCE, tag_seek, MPI_COMM_WORLD);
    
    actual_tsk++;

  }
  
  
  //MPI_Recv & MPI_Send not more work
  int numprocs;
  int actual_proc = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);  
  while (actual_proc < numprocs - 2) {
    
    MPI_Recv(&petition, 1, MPI_INT, MPI_ANY_SOURCE, tag_mng, MPI_COMM_WORLD, &status);

    s_positions[0] = -1;
    s_positions[1] = -1;

    MPI_Send(s_positions, 2*sizeof(size_t), MPI_CHAR, status.MPI_SOURCE, tag_seek, MPI_COMM_WORLD);

    actual_proc++;
  }
  

  printf("MANAGER FINISH!\n");

}


/*****************************************************************/
/*****************************************************************/
// N MPI PROCESS
// N - 1 : MPI WRITER PRCOESS
// N - 2 : MPI MANAGER PRCOESS
/*****************************************************************/
/*****************************************************************/

void rna_aligner_mpi(options_t *options, int argc, char *argv[]) {
  //MPI Process Files  in the nodes
  FILE* fd_out;

  int tag = MPI_TAG;
  int tag_mng  = MPI_TAG_MNG;
  int tag_seek = MPI_TAG_SEEK;

  MPI_Status status;
  int numprocs;
  int rank;
  int len;
  n_tasks_t ntasks;
  int th_steal = 0;
  int worker   = 1;
  int total_files;
  MPI_data_t MPI_data;

  size_t reads_wf1, reads_wf2, reads_wf3;

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;

  id_send = 0;
  
  mpi_tot_send = 0;
  t_mpi_send = 0;
  reads_wf = 0;
  reads_r = 0;
  
  pthread_mutex_init(&mutex_merge, NULL);
  
  MPI_data.block = 1;
  MPI_data.run = 1;
  pthread_mutex_init(&(MPI_data.queue_mutex), NULL);
  //Queue of tasks

  //printf("Init...\n");
  // ([0] Task for process this node) <--- [0]-[1]-[2]-[3] ---> ([3] Task for steal from thief node)
  MPI_data.queue_tasks = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);

  linked_list_t *queue_sa_tasks = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  linked_list_t *queue_hc_tasks = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);  
  
  int total_tasks = 0;
  int init_id;
  
  argc += 1;
  argv -= 1;
  
  MPI_Init(&argc,&argv);
  int provided;
  //MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  //MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  
  //if (provided != MPI_THREAD_MULTIPLE) {
  //printf("Sorry, this MPI implementation does not support multiple threads\n");
  //MPI_Abort(MPI_COMM_WORLD, 1);
  //}
  
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &name_len);

  //printf("[%i/%i] %s\n", rank, numprocs, processor_name);

  //================= S Y N C ==================//
  //unsigned char sync_data[numprocs];
  //memset(sync_data, 0, numprocs);
  //--------------------------------------------//
  
  //============== T I M I N G =================//
  double t_genomes, t_split, t_process, t_merge;               //t_file, t_file_write, t_total_file, t_loop;
  struct timeval start, stop, start_split, stop_split;         //start_file, stop_file, start_loop, stop_loop;

  double t_total = 0, t_w1 = 0, t_w2 = 0, t_w3 = 0, t_m1 = 0, t_m2 = 0, t_m3 = 0;

  struct timeval start_total, stop_total, start_w, stop_w;
  
  t_genomes = 0;
  t_split = 0;
  t_process = 0;
  t_merge = 0;

  double t_file;
  struct timeval start_file, stop_file;

  t_file = 0;
  //-------------------------------------------//

  
  //========================= F I L E S    P A T H S ==========================//
  //End fill 
  int path_length = strlen(options->output_name);
  int prefix_length = 0;
  
  if (options->prefix_name) {
    prefix_length = strlen(options->prefix_name);
  }
  
  char *reads_results = (char *)calloc((60 + prefix_length), sizeof(char));
  char *exact_junctions = (char *)calloc((60 + prefix_length), sizeof(char));
  char *output_filename = (char *)calloc((path_length + prefix_length + 60), sizeof(char));
  char *exact_filename = (char *)calloc((path_length + prefix_length + 60), sizeof(char));
  //struct timeval time_genome_s, time_genome_e;
  //double time_genome = 0.0;
  metaexons_t *metaexons;
  genome_t *genome;  
  bwt_index_t *bwt_index = NULL;
  sa_index3_t *sa_index = NULL;
  int num_chromosomes;
  // load index for dna/rna or for bisulfite case

  char filename_tab[strlen(options->bwt_dirname) + 1024];  
  sprintf(filename_tab, "%s/params.info", options->bwt_dirname);
  FILE *fd = fopen(filename_tab, "r");

  if (fd) { 
    options->fast_mode = 1; 
    fclose(fd);
  } else { 
    options->fast_mode = 0; 
  }

  
  if (options->fast_mode && !options->set_cal) {
    options->min_cal_size = 30;
  }

  if (options->fast_mode && options->set_cal) {
    if (options->min_cal_size < 30) {
      options->min_cal_size = 30;
    } else if(options->min_cal_size > 100)  {
      options->min_cal_size = 100;
    }
  }

  options->bam_format = 0;

  if (options->prefix_name) {
    strcat(reads_results, "/");
    strcat(reads_results, options->prefix_name);
    strcat(reads_results, "_");
    strcat(reads_results, OUTPUT_FILENAME);
    strcat(reads_results, ".sam");  
    strcat(exact_junctions, "/");
    strcat(exact_junctions, options->prefix_name);
    strcat(exact_junctions, "_exact_junctions.bed");
  } else {
    strcat(reads_results, "/");    
    strcat(reads_results, OUTPUT_FILENAME);
    strcat(reads_results, ".sam");      
    strcat(exact_junctions, "/exact_junctions.bed");
  } 

  strcat(output_filename, options->output_name);
  strcat(output_filename, reads_results);
  free(reads_results);

  strcat(exact_filename, options->output_name);
  strcat(exact_filename, exact_junctions);
  free(exact_junctions);

  char *file1 = options->in_filename;

  if (rank < numprocs - 2) {    
    printf("[%i] Load genome...\n", rank);
    if (!options->fast_mode) {
      //===================================================================//
      //                L O A D    S T R U C T U R E S                     //
      //===================================================================//
      //////////////// LOAD BWT INDEX //////////////////////    
      bwt_index = bwt_index_new(options->bwt_dirname, false);
      ////////////////     GENOME     //////////////////////    
      genome = genome_new("dna_compression.bin", options->bwt_dirname, BWT_MODE);  
      //===================================================================//
      //                                                                   //
      //===================================================================//
    } else {
      sa_index3_parallel_genome_new(options->bwt_dirname, options->num_cpu_threads, &sa_index, &genome);
    }
    printf("[%i] Load done!\n", rank);
    
    num_chromosomes = genome->num_chromosomes;
    
    // Metaexons structure
    metaexons = metaexons_new(genome->num_chromosomes, 
			      genome->chr_size);    
  }
  
  MPI_Barrier(MPI_COMM_WORLD);



  if (rank == 0) {
    //TODO: MOVE CORRECT POSITION
    create_directory(TMP_PATH_BUFFERS);  

    FILE* fd_p = fopen(output_filename, "w");
    write_sam_header_BWT(genome, (FILE *) fd_p);
    fclose(fd_p);
  }
   

  //MPI_Bcast(read_positions, sizeof(size_t)*numprocs, MPI_CHAR, 0, MPI_COMM_WORLD);

  if (rank == 0) { printf("START PROCESS...\n"); }
  
  int type_rec = 0;
  
  //printf("==========================================\n");
  
  bwt_optarg_t *bwt_optarg;
  cal_optarg_t *cal_optarg;
  pair_mng_t *pair_mng;
  report_optarg_t *report_optarg;
  fastq_batch_reader_input_t reader_input;
  linked_list_t *alignments_list;
  sw_optarg_t sw_optarg;  
  bwt_server_input_t bwt_input;
  region_seeker_input_t region_input;
  cal_seeker_input_t cal_input;
  pair_server_input_t pair_input;
  batch_writer_input_t writer_input;
  sw_server_input_t sw_input;
  avls_list_t* avls_list;
  sa_wf_batch_t *sa_wf_batch;
			    
  workflow_t *wf;    

  workflow_SA_t *sa_wf;
  wf_input_t *wf_input;
  sa_wf_input_t *sa_wf_input;
  
  FILE *f_sa, *f_hc;
  batch_t *batch;
  batch_buffer_t *data_out;

  sa_rna_input_t sa_rna;

  size_t msj_1, msj_2, msj_3;
  size_t mmeta_1, mmeta_2, mmeta_3;
  
  if (rank < numprocs - 2) {
    //==== T A S K S    I N S E R T S ====//
    //int item = init_id;
    //for (int it = 0; it < total_tasks; it++) {
    //char file_name[1024];
    //sprintf(file_name, "%i.tmp", item);
    //linked_list_insert_first(strdup(file_name), MPI_data.queue_tasks);
    //linked_list_insert_first(strdup(file_name), queue_sa_tasks);
    //linked_list_insert_first(strdup(file_name), queue_hc_tasks);
    //item++;
    //}
    
    //=================== I N P U T    I N I T I A L I Z A T I O N S =========================//
    //                                                                                        //
    //========================================================================================// 
    //BWT parameters
    extern int min_intron, max_intron;
    min_intron = options->min_intron_length;
    max_intron = options->max_intron_length;
  
    bwt_optarg = bwt_optarg_new(1, 0,
				options->filter_read_mappings, 
				options->filter_seed_mappings);  
    // CAL parameters
    //printf("%i\n", options->min_cal_size);
    cal_optarg = cal_optarg_new(options->min_cal_size, 
				options->seeds_max_distance, 
				options->num_seeds, 
				options->min_num_seeds_in_cal,
				options->seed_size, 
				options->min_seed_size, 
				options->cal_seeker_errors, 
				options->max_intron_length, 
				options->min_intron_length);

    // paired mode parameters
    pair_mng = pair_mng_new(options->pair_mode, options->pair_min_distance, 
			    options->pair_max_distance, options->report_only_paired);

    // report parameters
    report_optarg = report_optarg_new(options->report_all,
				      options->report_n_best,
				      options->report_n_hits, 
				      options->report_only_paired,
				      options->report_best);


    linked_list_t *buffer    = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);
    linked_list_t *buffer_hc = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);
    
    fastq_batch_reader_input_init(options->in_filename, options->in_filename2, 
				  options->pair_mode, options->batch_size, 
				  NULL, options->gzip, &reader_input);  

    alignments_list = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);
    linked_list_set_flag(options->pair_mode, alignments_list);

    sw_optarg_init(options->gap_open, options->gap_extend, 
		   options->match, options->mismatch, &sw_optarg);

    bwt_server_input_init(NULL,  0,  bwt_optarg, 
			  bwt_index, NULL,  0, 
			  NULL, metaexons, &sw_optarg, 
			  genome, &bwt_input);
      
    region_seeker_input_init(NULL, cal_optarg, 
			     bwt_optarg, bwt_index, NULL, 
			     0, 0, options->min_seed_padding_left, 
			     options->min_seed_padding_right, genome, metaexons, &region_input);
  
    cal_seeker_input_init(NULL, cal_optarg, NULL, 0, NULL, NULL, genome, 
			  bwt_optarg, bwt_index, metaexons, &cal_input);
  
    int pair_mode = pair_mng->pair_mode;

    avls_list = avls_list_new(num_chromosomes);

    sw_server_input_init(NULL, NULL, 0,  options->match,  
			 options->mismatch,  options->gap_open, options->gap_extend,  
			 options->min_score,  options->flank_length, genome,  
			 options->max_intron_length, options->min_intron_length,  
			 options->seeds_max_distance,  bwt_optarg, avls_list, 
			 cal_optarg, bwt_index, metaexons, buffer, buffer_hc, 
			 NULL, NULL, pair_mode, &sw_input);
    
    pair_server_input_init(pair_mng, report_optarg, NULL, NULL, NULL, &pair_input);

    batch_writer_input_init(output_filename,
			    exact_filename,
			    NULL,
			    //extend_filename, 
			    alignments_list, 
			    genome, 
			    &writer_input);

    writer_input.bam_format = options->bam_format;
    writer_input.bam_file = (bam_file_t *) fopen(output_filename, "w"); 

    //clear_batch_buffer(batch_buffer_t *bb)

    data_out = new_batch_buffer();
        
    //Parse input for more than one file

    fastq_batch_reader_input_init(file1, NULL,
				  options->pair_mode, options->batch_size, 
				  NULL, options->gzip, &reader_input);  

    fflush(stdout);

    if (options->fast_mode) {
      //sw_optarg_t sw_optarg;
      //sw_optarg_init(options->gap_open, options->gap_extend, 
      //	     options->match, options->mismatch, &sw_optarg);      
      sa_rna.cal_optarg   = cal_optarg;
      sa_rna.genome    = genome;
      sa_rna.avls_list = avls_list;
      sa_rna.metaexons = metaexons;
      sa_rna.sw_optarg = &sw_optarg;
      sa_rna.file1 = f_sa;
      sa_rna.file2 = f_hc;
      sa_rna.pair_input = &pair_input;
      sa_rna.min_score = options->min_score;
      sa_rna.max_alig = options->filter_read_mappings;
      
      sa_wf_batch = sa_wf_batch_new(NULL, (void *)sa_index, &writer_input, NULL, &sa_rna);
      sa_wf_batch->buffer_mpi = data_out;
      sa_wf_input = sa_wf_input_new(options->bam_format, &reader_input, sa_wf_batch);
      
      // create and initialize workflow
      sa_wf = workflow_SA_new(); 
      workflow_stage_function_SA_t stage_functions[] = {sa_rna_mapper};
      char *stage_labels[] = {"SA mapper"};
      workflow_set_stages_SA(1, stage_functions, stage_labels, sa_wf);      
      // optional producer and consumer functions      
      //workflow_set_producer_SA((workflow_producer_function_SA_t *)sa_fq_reader_rna, "FastQ reader", sa_wf);
      workflow_set_producer_SA((workflow_producer_function_SA_t *)SA_MPI_producer, "FastQ reader", sa_wf);      
      workflow_set_consumer_SA((workflow_consumer_function_SA_t *)SA_MPI_consumer, "SAM writer", sa_wf);      
    } else {    
      //init_data_output(&data_out, READS_SPLIT);
      batch = batch_new(&bwt_input, &region_input, &cal_input, 
			&pair_input, NULL, &sw_input, &writer_input, RNA_MODE, NULL, data_out);
      
      wf_input = wf_input_new(&reader_input, batch);  
      //create and initialize workflow
      wf = workflow_new();
      //workflow_stage_function_t stage_functions[] = {bwt_stage, cal_stage, 
      //				     sw_stage, post_pair_stage, MPI_output};  
      //char *stage_labels[] = {"BWT", "CAL", "SW", "POST PAIR", "MPI OUTPUT"};
      //workflow_set_stages(5, (workflow_stage_function_t *)&stage_functions, stage_labels, wf);

      workflow_stage_function_t stage_functions[] = {bwt_stage, cal_stage, 
						     sw_stage, post_pair_stage};
      char *stage_labels[] = {"BWT", "CAL", "SW", "POST PAIR"};
      workflow_set_stages(4, (workflow_stage_function_t *)&stage_functions, stage_labels, wf);

      //optional producer and consumer functions
      //workflow_set_producer((workflow_producer_function_t *)fastq_reader, "FastQ reader", wf);
      workflow_set_producer((workflow_producer_function_t *)BWT_MPI_producer, "FastQ reader", wf); 
      workflow_set_consumer((workflow_consumer_function_t *)MPI_output, "SAM writer", wf);
      //workflow_set_consumer((workflow_consumer_function_t *)sam_writer, "SAM writer", wf);
    }
  }
  
  //MPI_Barrier(MPI_COMM_WORLD);
  
  //if (rank == 0) {
  //printf("\t START MAPPING...\n");
  //start_timer(start_w);
  //}

  double tot_t_file = 0;
  double tot_t_send = 0;

  mpi_tot_reads = 0;
  r_tot_reads = 0;

  struct timeval wx_start, wx_stop;
  double wx_time;
  double t_merge_1;
  
  if (rank == 0) {
    start_timer(wx_start);
  }
  

  
  //create new communicator
  //--------------------------------------------------------------  
  MPI_Comm new_comm; 
  MPI_Group orig_group, new_group;
  MPI_Comm_group(MPI_COMM_WORLD, &orig_group); 
  
  /* Divide tasks into two distinct groups based upon rank */ 
  //if (rank != numprocs -1) {
  int ranks_group[numprocs - 2];
  for (int i = 0; i < numprocs - 2; i++) {
    ranks_group[i] = i;
  }
  MPI_Group_incl(orig_group, numprocs - 2, ranks_group, &new_group);
  
  /* Create new communicator and then perform collective communications */
  MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm); 
  
  //--------------------------------------------------------------

  //-> CALL MANAGER
  if (rank == numprocs - 2) {
    manager_function(file1);
  } else if (rank < numprocs - 2) {     

    //===================================================================//
    //**********************      NEW MPI      **************************//
    //===================================================================//

    size_t s_positions[2];
    int petition = 1;
    MPI_Status status;
    while (1) {
      
      MPI_Send(&petition, 1, MPI_INT, numprocs - 2, tag_mng, MPI_COMM_WORLD);   
      MPI_Recv(s_positions, 2*sizeof(size_t), MPI_CHAR, numprocs - 2, tag_seek, MPI_COMM_WORLD, &status);
      
      size_t start_seek = s_positions[0];
      size_t end_seek   = s_positions[1];

      printf("[%i]SEEK: %lld - %lld\n", rank, start_seek, end_seek);
      sleep(1);
      
      if ((start_seek == -1) && (end_seek == -1)) { break; }

      
    }

    printf("[%i]WORKER FINISH!\n", rank);

    //===================================================================//
    //*******************************************************************//
    //===================================================================//

    /*
    size_t num_cpus = get_optimal_cpu_num_threads();    
    //==== C R E A T E    P O L L I N G    T H R E A D ====//    
    //=============================//
    //= W O R K E R - T H R E A D =//
    //=============================//
    char *file;
    //int victim;
    //work_package_t work;

    double t_task = 0, t_send = 0;
    struct timeval start_task, stop_task, start_send, stop_send;
    int w1_tasks = 0;
    

    //==============-----> START loop

    //TODO: MPI Send -> 
    //TODO: MPI Recv <-

    //====       W O R K E R    P R O C E S S    F I L E       ====//
    //=============================================================//
    
    char file_path[1024];
    char buffer_hc_path[1024];
    
    //sprintf(file_path, "%s/%s", TMP_PATH_FILES, file);
    sprintf(buffer_hc_path, "%s/buffer_hc.%i", TMP_PATH_BUFFERS, rank);

    start_seek = read_positions[rank];
    end_seek   = read_positions[rank + 1];
    //pos_seek   = start_seek;
    
    reader_input.fq_file1 = fastq_fopen(file1);
    fseek(reader_input.fq_file1->fd, start_seek, SEEK_SET);
    
    f_hc = fopen(buffer_hc_path, "w+b");
    if (f_hc == NULL) {
      LOG_FATAL("Error opening file 'buffer_hc.tmp' \n");
    }
    
    t_file = 0;
    if (options->fast_mode) {
      sa_rna.file1 = f_hc;
      struct timeval start_SA, stop_SA;
      double t_SA = 0;
      start_timer(start_SA);
      workflow_run_with_SA(options->num_cpu_threads, sa_wf_input, sa_wf);
      stop_timer(start_SA, stop_SA, t_SA);

      //printf("[%i]Process SA file with %i threads in %0.2f, total(s) reads = %i\n", rank, options->num_cpu_threads, t_SA / 1000000, r_tot_reads);
    
      sa_wf->completed_producer = 0;
    } else {
      struct timeval start_BWT, stop_BWT;
      double t_BWT = 0;
      
      char buffer_sa_path[1024];
      sprintf(buffer_sa_path, "%s/buffer_sa.%i", TMP_PATH_BUFFERS, rank);      

      f_sa = fopen(buffer_sa_path, "w+b");
      if (f_sa == NULL) {
	LOG_FATAL("Error opening file 'buffer_hc.tmp' \n");
      }

      sw_input.f_sa = f_sa;
      sw_input.f_hc = f_hc;

      //printf("Start BWT\n");
      start_timer(start_BWT);
      //for (int k = 0; k < 100; k++) {
      //size_t *kaka = (size_t *)malloc(sizeof(size_t)*1000000);
      //free(kaka);
      //}
      workflow_run_with(options->num_cpu_threads, wf_input, wf);
      workflow_display_timing(wf);
      
      stop_timer(start_BWT, stop_BWT, t_BWT);
      wf->completed_producer = 0;
      //printf("[%i]Process BWT file with %i threads in %0.2f, total(s) reads = %i\n", rank, options->num_cpu_threads, t_BWT / 1000000, reads_r);
      fclose(f_sa); 
    }

    fclose(f_hc); 

    //==============-----> END loop
    
    fastq_fclose(reader_input.fq_file1);    
        
    batch_out_t *b = array_list_get(1, data_out->buffer);    
    //Send Batches to Write      
    MPI_Send(b->data, MAX_BUFFER, MPI_CHAR, numprocs - 1, 1, MPI_COMM_WORLD);
    b->len = 0;




    

    //==========================================================================//
    //= M E R G E    M E T A E X O N - A V L     T O    2nd    W O R K F L O W =//
    //==========================================================================//
    printf("[%i]Start merge...\n", rank);
    MPI_Barrier(new_comm);


    struct timeval start_BWT, stop_BWT;
    start_timer(start_BWT);
    
    unsigned long num_sj;        
    int numprocs_merge = numprocs - 1;
    int dec_s = 1;
    int dec_r = 2;
    int send, recv;
    int dec = 2;
    
    while (1) {
      send = numprocs_merge - dec_s;
      recv = numprocs_merge - dec_r;
      if (send > 0 && recv < 0)
	recv = 0;
      if (send <= 0 && recv <= 0)
	break;
      
      for (int r = numprocs_merge - 1; r >= 0; r--) {
	if (rank == send) {
	  //SEND AVL
	  num_sj = get_total_items(num_chromosomes, avls_list);
	  unsigned long  *list_sj = (unsigned long *)malloc(sizeof(unsigned long)*num_sj*5); 	  
	  MPI_avl_package(num_chromosomes, avls_list, num_sj, list_sj);	  
	  
	  
	  //SEND METAEXON
	  unsigned long num_meta = 0;
	  unsigned long num_left_breaks = 0, num_right_breaks = 0;
	  for (int c = 0; c < metaexons->num_chromosomes; c++) {
	    linked_list_t *metaexon_list = metaexons->metaexons_list[c];
	    linked_list_item_t *item;
	    for (item = metaexon_list->first; item != NULL; item = item->next) {
	      metaexon_t *metaexon = item->item;
	      num_left_breaks  += array_list_size(metaexon->left_breaks);
	      num_right_breaks += array_list_size(metaexon->right_breaks);
	      num_meta++;	      
	    }
	  }
	  
	  //printf("Pack meta %i, %i, %i...\n", num_meta, num_left_breaks, num_right_breaks);
	  unsigned long *list_meta = (unsigned long *)malloc(sizeof(unsigned long)*num_meta*5);
	  unsigned long *list_left_breaks  = (unsigned long *)malloc(sizeof(unsigned long)*num_left_breaks*2);
	  unsigned long *list_right_breaks = (unsigned long *)malloc(sizeof(unsigned long)*num_right_breaks*2);
	  
	  MPI_metaexon_package(metaexons, num_meta, list_meta, 
			       num_left_breaks, list_left_breaks, 
			       num_right_breaks, list_right_breaks);
	  
	  //printf("RECV-RANK0: %i, %i, %i, %i, %i\n", list_meta[num_meta], list_meta[num_meta + 1], list_meta[num_meta + 2], list_meta[num_meta + 3], list_meta[num_meta + 4]);

	  MPI_Send(&num_sj, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);	  
	  MPI_Send(&num_meta, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(&num_left_breaks, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(&num_right_breaks, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  //printf("[%i]SEND VALUES : %lu, %lu, %lu, %lu\n", rank, num_sj, num_meta, num_left_breaks, num_right_breaks);	  

	  MPI_Send(list_sj, num_sj*5, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(list_meta, num_meta*5, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(list_left_breaks, num_left_breaks*2, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(list_right_breaks, num_right_breaks*2, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);

	  //printf("[%i]SEND LISTS : %lu, %lu, %lu, %lu\n", rank, list_sj[num_sj], list_meta[num_meta], list_sj[num_sj + 1], list_meta[num_meta + 1]);
	  //printf("SEND: %i, %i, %i, %i, %i\n", list_right_breaks[0], list_right_breaks[1], list_right_breaks[2], list_right_breaks[3], list_right_breaks[4]);
	
	  free(list_sj);
	  free(list_meta);
	  free(list_left_breaks);
	  free(list_right_breaks);
	
	}
	
	if (rank == recv) {
	  //Recv AVLs
	  unsigned long recv_num_sj;	    
	  unsigned long recv_num_meta;
	  unsigned long recv_num_left_breaks;
	  unsigned long recv_num_right_breaks;
	
	  //printf("RECV: %i, %i, %i, %i, %i\n", recv_list_meta[recv_num_meta], recv_list_meta[recv_num_meta + 1], recv_list_meta[recv_num_meta + 2], recv_list_meta[recv_num_meta + 3], recv_list_meta[recv_num_meta + 4]);	  
	
	  MPI_Recv(&recv_num_sj, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);	
	  MPI_Recv(&recv_num_meta, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  MPI_Recv(&recv_num_left_breaks, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  MPI_Recv(&recv_num_right_breaks, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  //printf("[%i]RECV VALUES : %lu, %lu, %lu, %lu\n", rank, recv_num_sj, recv_num_meta, recv_num_left_breaks,recv_num_right_breaks);	  

	  unsigned long *recv_list_sj = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_sj*5);
	  unsigned long *recv_list_meta = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_meta*5);
	  unsigned long *recv_list_left_breaks = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_left_breaks*2);
	  unsigned long *recv_list_right_breaks = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_right_breaks*2);

	  MPI_Recv(recv_list_sj, recv_num_sj*5, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);	  	
	  MPI_Recv(recv_list_meta, recv_num_meta*5, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  MPI_Recv(recv_list_left_breaks, recv_num_left_breaks*2, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  MPI_Recv(recv_list_right_breaks, recv_num_right_breaks*2, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);


	  //Merge AVLs
	  avl_node_t *node_avl_start, *node_avl_end;
	  MPI_splice_t MPI_sj;
	  for (size_t sj = 0; sj < recv_num_sj; sj++) {
	    MPI_sj.start          = recv_list_sj[sj];
	    MPI_sj.end            = recv_list_sj[recv_num_sj + sj];
	    MPI_sj.reads_number   = recv_list_sj[recv_num_sj*2 + sj];
	    MPI_sj.strand         = recv_list_sj[recv_num_sj*3 + sj];
	    MPI_sj.chr            = recv_list_sj[recv_num_sj*4 + sj];
	  
	    //printf("\t[%i] %i/%i Merge (%i) %i: %lu - %lu\n", rank, sj, recv_num_sj, MPI_sj.strand, MPI_sj.chr, MPI_sj.start, MPI_sj.end);
	    allocate_start_node(MPI_sj.chr - 1,
				MPI_sj.strand,
				MPI_sj.start,
				MPI_sj.end,
				MPI_sj.start,
				MPI_sj.end,
				FROM_READ,
				1,
				"",
				&node_avl_start,
				&node_avl_end,
				avls_list);
	    //Refresh MPI_sj.reads_number in nodes
	    //end_data_t *end_data = node_avl_start->data;	    
	  }
	  //Merge Metaexon
	  //sleep(1);
	
	  //printf("[%i]Merge %i META:\n", rank, recv_num_meta);
	  size_t l_pos = 0, r_pos = 0;
	  MPI_metaexon_t mpi_meta;
	  for (size_t m = 0; m < recv_num_meta; m++) {
	    mpi_meta.chromosome = recv_list_meta[m];
	    mpi_meta.start      = recv_list_meta[recv_num_meta   + m];
	    mpi_meta.end        = recv_list_meta[recv_num_meta*2 + m];
	    mpi_meta.n_starts   = recv_list_meta[recv_num_meta*3 + m];
	    mpi_meta.n_ends     = recv_list_meta[recv_num_meta*4 + m];

	    if (!mpi_meta.n_starts && !mpi_meta.n_ends) {
	      //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);	      
	      metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			      METAEXON_LEFT_END, NULL, metaexons);
	    }
	    
	    
	    MPI_breaks_t bk;
	    for (size_t s = 0; s < mpi_meta.n_starts; s++) {
	      bk.pos    = recv_list_right_breaks[r_pos];
	      bk.strand = recv_list_right_breaks[recv_num_right_breaks + r_pos];
	      r_pos++;
	      //MPI_breaks_t bk = MPI_breaks[se_pos++];
	      //printf("\tM-STARTS-(%i)%lu\n", bk.strand, bk.pos);
	      cp_avltree *avl = avls_list->avls[bk.strand][mpi_meta.chromosome].avl;
	      avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	      assert(avl_node);
	      //if (!avl_node) {
	      //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
	      //printf("%i:%i:%i\n", m, r_pos - 1, bk.pos);
	      //exit(-1);
	      //}
	      metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			      METAEXON_RIGHT_END, avl_node, metaexons);
	    }
	    for (size_t s = 0; s < mpi_meta.n_ends; s++) {
	      bk.pos    = recv_list_left_breaks[l_pos];
	      bk.strand = recv_list_left_breaks[recv_num_left_breaks + l_pos];
	      l_pos++;
	      //MPI_breaks_t bk = MPI_breaks[se_pos++];
	      //printf("\t[%i]M-ENDS-(%i)%lu\n", rank, bk.strand, bk.pos);
	      cp_avltree *avl = avls_list->ends_avls[bk.strand][mpi_meta.chromosome].avl;
	      avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	      assert(avl_node);
	      metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			      METAEXON_LEFT_END, avl_node, metaexons);
	    }
	  }
	
	  free(recv_list_sj);
	  free(recv_list_meta);
	  free(recv_list_left_breaks);
	  free(recv_list_right_breaks);

	}
	send -= dec;
	recv -= dec;
	if (send <= 0) { send = -1; }
      }
      dec   *= 2;
      dec_s *= 2;
      dec_r *= 2;
    }
  
    //printf("[%i]MERGE NODO 0 WORKFLOW 1 END...BROADCAST TO OTHER NODES\n", rank);
  
    //int num_sj;
    unsigned long recv_num_sj, recv_n_metaexons, recv_n_starts_ends;
  
    unsigned long *list_sj;
    unsigned long *list_meta;
    unsigned long *list_left_breaks;
    unsigned long *list_right_breaks;
    
    unsigned long num_meta = 0;
    unsigned long num_left_breaks = 0, num_right_breaks = 0;    
    
    if (rank == 0) {
      //PACK AVL
      num_sj = get_total_items(num_chromosomes, avls_list);
      list_sj = (unsigned long *)malloc(sizeof(unsigned long)*num_sj*5); 	  
      MPI_avl_package(num_chromosomes, avls_list, num_sj, list_sj);	  
      
      //PACK METAEXON
      for (int c = 0; c < metaexons->num_chromosomes; c++) {
	linked_list_t *metaexon_list = metaexons->metaexons_list[c];
	linked_list_item_t *item;
	for (item = metaexon_list->first; item != NULL; item = item->next) {
	  metaexon_t *metaexon = item->item;
	  num_left_breaks  += array_list_size(metaexon->left_breaks);
	  num_right_breaks += array_list_size(metaexon->right_breaks);
	  num_meta++;	      
	}
      }
      
      list_meta = (unsigned long *)malloc(sizeof(unsigned long)*num_meta*5);
      list_left_breaks  = (unsigned long *)malloc(sizeof(unsigned long)*num_left_breaks*2);
      list_right_breaks = (unsigned long *)malloc(sizeof(unsigned long)*num_right_breaks*2);

      MPI_metaexon_package(metaexons, num_meta, list_meta, 
			   num_left_breaks, list_left_breaks, 
			   num_right_breaks, list_right_breaks);
    }
    
    MPI_Bcast(&num_sj, 1, MPI_UNSIGNED_LONG, 0, new_comm);
    MPI_Bcast(&num_meta, 1, MPI_UNSIGNED_LONG, 0, new_comm);
    MPI_Bcast(&num_left_breaks, 1, MPI_UNSIGNED_LONG, 0, new_comm);
    MPI_Bcast(&num_right_breaks, 1, MPI_UNSIGNED_LONG, 0, new_comm);

    size_t xn = num_sj*5 + num_meta*5 + num_left_breaks*2 + num_right_breaks*2;
    unsigned long *pack_send = (unsigned long *)malloc(sizeof(unsigned long)*xn);

    if (rank == 0) {
      memcpy(pack_send, list_sj, sizeof(unsigned long)*num_sj*5);
      memcpy(&pack_send[num_sj*5], list_meta, sizeof(unsigned long)*num_meta*5);
      memcpy(&pack_send[num_sj*5 + num_meta*5], list_left_breaks, sizeof(unsigned long)*num_left_breaks*2);
      memcpy(&pack_send[num_sj*5 + num_meta*5 + num_left_breaks*2], list_right_breaks, sizeof(unsigned long)*num_right_breaks*2);

      free(list_sj);
      free(list_meta);
      free(list_left_breaks);
      free(list_right_breaks);
    }

    MPI_Bcast(pack_send, xn, MPI_UNSIGNED_LONG, 0, new_comm);

    list_sj           = pack_send;
    list_meta         = &pack_send[num_sj*5];
    list_left_breaks  = &pack_send[num_sj*5 + num_meta*5];
    list_right_breaks = &pack_send[num_sj*5 + num_meta*5 + num_left_breaks*2];

    //Merge AVLs
    avl_node_t *node_avl_start, *node_avl_end;
    MPI_splice_t MPI_sj;

    for (size_t sj = 0; sj < num_sj; sj++) {
      MPI_sj.start          = list_sj[sj];
      MPI_sj.end            = list_sj[num_sj + sj];
      MPI_sj.reads_number   = list_sj[num_sj*2 + sj];
      MPI_sj.strand         = list_sj[num_sj*3 + sj];
      MPI_sj.chr            = list_sj[num_sj*4 + sj];

      //printf("[%i]Merge 2 (%i) %i: %lu - %lu\n", rank, MPI_sj.strand, MPI_sj.chr, MPI_sj.start, MPI_sj.end);
      allocate_start_node(MPI_sj.chr - 1,
			  MPI_sj.strand,
			  MPI_sj.start,
			  MPI_sj.end,
			  MPI_sj.start,
			  MPI_sj.end,
			  FROM_READ,
			  1,
			  "",
			  &node_avl_start,
			  &node_avl_end,
			  avls_list);
      //Refresh MPI_sj.reads_number in nodes
      end_data_t *end_data = node_avl_start->data;	    
    }
    
    //printf("[%i]Merge Meta %i\n", rank, num_meta);
    
    //Merge Metaexon
    size_t se_pos = 0;
    size_t r_pos = 0, l_pos = 0;
    MPI_metaexon_t mpi_meta;
    
    for (size_t m = 0; m < num_meta; m++) {
      mpi_meta.chromosome = list_meta[m];
      mpi_meta.start      = list_meta[num_meta   + m];
      mpi_meta.end        = list_meta[num_meta*2 + m];
      mpi_meta.n_starts   = list_meta[num_meta*3 + m];
      mpi_meta.n_ends     = list_meta[num_meta*4 + m];

      if (!mpi_meta.n_starts && !mpi_meta.n_ends) {
	//printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends)
	metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			METAEXON_LEFT_END, NULL, metaexons);
      } 
      
      //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
      MPI_breaks_t bk;
      for (int s = 0; s < mpi_meta.n_starts; s++) {
	bk.pos    = list_right_breaks[r_pos];
	bk.strand = list_right_breaks[num_right_breaks + r_pos];
	r_pos++;

	//printf("\tM-STARTS-(%i)%lu\n", bk.strand, bk.pos);
	cp_avltree *avl = avls_list->avls[bk.strand][mpi_meta.chromosome].avl;
	avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	//assert(avl_node);
	if (!avl_node) {
	  printf("ERROR: (%i|%i)%lu\n", m, s, bk.pos);
	  exit(-1);
	}
	metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			METAEXON_RIGHT_END, avl_node, metaexons);
      }
      for (int s = 0; s < mpi_meta.n_ends; s++) {
	bk.pos    = list_left_breaks[l_pos];
	bk.strand = list_left_breaks[num_left_breaks + l_pos];
	l_pos++;

	//printf("\t[%i]M-ENDS-(%i)%lu\n", rank, bk.strand, bk.pos);
	cp_avltree *avl = avls_list->ends_avls[bk.strand][mpi_meta.chromosome].avl;
	avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	assert(avl_node);
	metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			METAEXON_LEFT_END, avl_node, metaexons);
      }
    }
  
    //printf("Frees merge\n");
    free(pack_send);
  
    //====                E N D    M E R G E                ====//
    //==========================================================//
    MPI_Barrier(new_comm);

    t_merge_1 = 0;
    stop_timer(start_BWT, stop_BWT, t_merge_1);
    
    //if (rank == 0) {
    //printf("Merge 1 time: %0.2f\n", t_BWT / 1000000);
    //}
    */
    //End Writer
    char *end = (char *)malloc(sizeof(char)*MAX_BUFFER);
    strcpy(end, "END\0");    
    
    MPI_Send(end, MAX_BUFFER, MPI_CHAR, numprocs - 1, 1, MPI_COMM_WORLD);
    free(end);      
    
  } else {
    
    double t_write = 0;
    struct timeval start_task, stop_task;
    
    //printf("[%i]MPI THREAD WRITER START\n", rank);
    
    //=============================//
    //= W R I T E R - T H R E A D =//
    //=============================//
    fd_out = fopen(output_filename, "a");
    
    if (type_rec == 0) {
      char *recv_buff = (char *)malloc(sizeof(char)*MAX_BUFFER);
      int numprocs_end = 0;
      while (1) {
	MPI_Recv(recv_buff, MAX_BUFFER, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
	if (strcmp("END", recv_buff) == 0) {
	  numprocs_end++;
	  if (numprocs_end == numprocs - 2) {
	    break;
	  }
	} else {
	  fwrite(recv_buff, strlen(recv_buff), sizeof(char), fd_out);
	}    
      }

    printf("WRITER FINISH!\n");

    } else {
      
      linked_list_t *writer_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
      pthread_cond_t writer_cond;
      pthread_mutex_t writer_mutex;
      int end_recv = 0;
      pthread_mutex_init(&writer_mutex, NULL);
      
      #pragma omp parallel sections num_threads(2)
      {	
          #pragma omp section
          {
	    char *recv_buff;
	    int numprocs_end = 0;	    

	    while (1) {
	      recv_buff = (char *)malloc(sizeof(char)*MAX_BUFFER);
	      MPI_Recv(recv_buff, MAX_BUFFER, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
	      if (strcmp("END", recv_buff) == 0) {
		numprocs_end++;
		if (numprocs_end == numprocs - 1) {
		  pthread_mutex_lock(&writer_mutex);
		  //printf("Insert item END\n");
		  linked_list_insert(recv_buff, writer_list);
		  end_recv = 1;
		  pthread_cond_signal(&writer_cond);
		  pthread_mutex_unlock(&writer_mutex);		  
		  break;
		}
		free(recv_buff);
	      } else {
		pthread_mutex_lock(&writer_mutex);
		//printf("Insert item recv\n");
		linked_list_insert(recv_buff, writer_list);
		pthread_cond_signal(&writer_cond);
		pthread_mutex_unlock(&writer_mutex);
	      }
	    }
	    //printf("Thread MPI Recv END\n");
	  }
	  #pragma omp section
          {
	    char *recv_buff;
	    while (1) {
	      pthread_mutex_lock(&writer_mutex);
	      while (!linked_list_size(writer_list) && !end_recv) {
		pthread_cond_wait(&writer_cond, &writer_mutex);
	      }
	      recv_buff = linked_list_remove_last(writer_list);
	      //printf("\tRemove item\n");
	      pthread_mutex_unlock(&writer_mutex);

	      if (strcmp("END", recv_buff) == 0) {
		//printf("\tItem END\n");
		free(recv_buff);
		break;
	      }
	      
	      fwrite(recv_buff, strlen(recv_buff), sizeof(char), fd_out);
	      free(recv_buff);
	    }
	    //printf("Thread MPI Writer END\n");
	  }
      }
    }    
    
  }

  //TODO: END WORK, DELETE THIS //
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();  

  return;

  /*
  if (rank == 0) {
    stop_timer(wx_start, wx_stop, wx_time);
    printf("============================WORKFLOW 1 TIME %0.2fs=========================\n", wx_time / 1000000);
    wx_time = 0;
    start_timer(wx_start);
  }
  */

  if (options->fast_mode) { goto fast_mode_skip; }

  double t_task = 0, t_send = 0;
  struct timeval start_task, stop_task, start_send, stop_send;
  double t_merge_2;
  
  if (rank != numprocs - 1) {
    //====            W2 - W O R K I N G               ====//
    //==== C R E A T E    P O L L I N G    T H R E A D ====//
    MPI_data.queue_tasks = queue_sa_tasks;
    wf_input_file_t *wf_input_file = wf_input_file_new(NULL, batch);   

    workflow_t *wf_last = workflow_new();
    workflow_stage_function_t stage_functions_last[] = {rna_last_stage, post_pair_stage};
    char *stage_labels_last[] = {"RNA LAST STAGE", "POST PAIR"};
    workflow_set_stages(2, (workflow_stage_function_t *)&stage_functions_last, stage_labels_last, wf_last);
    workflow_set_producer((workflow_producer_function_t *)file_reader, "Buffer reader", wf_last);
    //workflow_set_consumer((workflow_consumer_function_t *)MPI_consumer, "SAM writer", wf_last);
    workflow_set_consumer((workflow_consumer_function_t *)MPI_output, "SAM writer", wf_last);
    
    //size_t num_cpus = get_optimal_cpu_num_threads();
    //num_cpus -= 2;

    //=============================//
    //= W O R K E R - T H R E A D =//
    //=============================//
  
    char *file;
    int victim;
    work_package_t work;

    start_timer(start_file);
    //====       W O R K E R    P R O C E S S    F I L E       ====//
    //=============================================================//
    char file_path[1024];
    char buffer_sa_path[1024];
    char buffer_hc_path[1024];
    
    //sprintf(file_path, "%s/buffer_sa.%s", TMP_PATH_FILES, file);
    sprintf(buffer_sa_path, "%s/buffer_sa.%i", TMP_PATH_BUFFERS, rank);
    sprintf(buffer_hc_path, "%s/buffer_hc.%i", TMP_PATH_BUFFERS, rank);
    
    //printf("rank.%i process file -> %s\n", rank, file_path);
    
    //Run workflow
    f_sa = fopen(buffer_sa_path, "r+b");
    if (f_sa == NULL) {
      LOG_FATAL("Error opening file 'buffer_sa.tmp' \n");
    }
    sw_input.f_sa = f_sa;
    
    f_hc = fopen(buffer_hc_path, "a+b");
    if (f_hc == NULL) {
      LOG_FATAL("Error opening file 'buffer_hc.tmp' \n");
    }
    
    wf_input_file->file = f_sa;
    
    //printf("[%i]Process file %s\n", rank, buffer_sa_path);
    start_timer(start_task);
    sw_input.f_hc = f_hc;
    t_file = 0;
    //printf("===== STOP ====\n");
    workflow_run_with(options->num_cpu_threads, wf_input_file, wf_last);
    //workflow_display_timing(wf_last);
    
    //printf("===== START ====\n");
    //workflow_run_with(1, wf_input_file, wf_last);
    //printf("===== STOP ====\n");
    tot_t_file += t_file;
    stop_timer(start_task, stop_task, t_task);
    
    wf_last->completed_producer = 0;
    
    fclose(f_sa);
    fclose(f_hc);
    
    //==== W O R K E R    E N D    P R O C E S S    F I L E ====//
    //==========================================================//
    stop_timer(start_file, stop_file, t_file);
    //printf("\t[%i]Task '%s' completed in %0.2f\n", rank, buffer_sa_path, t_file / 1000000);
    t_file = 0;
    //}
    
    batch_out_t *b = array_list_get(1, data_out->buffer);    
    //printf("[%i]SEND  MPI\n", rank);
    MPI_Send(b->data, MAX_BUFFER, MPI_CHAR, numprocs - 1, 1, MPI_COMM_WORLD);
    b->len = 0;


    //==========================================================================//
    //= M E R G E    M E T A E X O N - A V L     T O    3nd    W O R K F L O W =//
    //==========================================================================//
    MPI_Barrier(new_comm);

    //////////////////////////////////////////////////////////////////
    /*
    {
      unsigned long num_sj = get_total_items(num_chromosomes, avls_list);
      //SEND METAEXON
      unsigned long num_meta = 0;
      unsigned long num_left_breaks = 0, num_right_breaks = 0;
      for (int c = 0; c < metaexons->num_chromosomes; c++) {
	linked_list_t *metaexon_list = metaexons->metaexons_list[c];
	linked_list_item_t *item;
	for (item = metaexon_list->first; item != NULL; item = item->next) {
	  metaexon_t *metaexon = item->item;
	  num_left_breaks  += array_list_size(metaexon->left_breaks);
	  num_right_breaks += array_list_size(metaexon->right_breaks);
	  num_meta++;	      
	}
      }
      printf("2.[%i]num_sj = %lu, num_meta = %lu\n", rank, num_sj, num_meta);
      msj_2 = num_sj;
      mmeta_2 = num_meta;
    }
    MPI_Barrier(new_comm);
    */
    //////////////////////////////////////////////////////////////////

    struct timeval start_BWT_2, stop_BWT_2;
    start_timer(start_BWT_2);

    unsigned long num_sj;        
    int numprocs_merge = numprocs - 1;
    int dec_s = 1;
    int dec_r = 2;
    int send, recv;
    int dec = 2;
    
    while (1) {
      send = numprocs_merge - dec_s;
      recv = numprocs_merge - dec_r;
      if (send > 0 && recv < 0)
	recv = 0;
      if (send <= 0 && recv <= 0)
	break;
      
      for (int r = numprocs_merge - 1; r >= 0; r--) {
	if (rank == send) {
	  //SEND AVL
	  num_sj = get_total_items(num_chromosomes, avls_list);
	  unsigned long  *list_sj = (unsigned long *)malloc(sizeof(unsigned long)*num_sj*5); 	  
	  MPI_avl_package(num_chromosomes, avls_list, num_sj, list_sj);	  


	  //SEND METAEXON
	  unsigned long num_meta = 0;
	  unsigned long num_left_breaks = 0, num_right_breaks = 0;
	  for (int c = 0; c < metaexons->num_chromosomes; c++) {
	    linked_list_t *metaexon_list = metaexons->metaexons_list[c];
	    linked_list_item_t *item;
	    for (item = metaexon_list->first; item != NULL; item = item->next) {
	      metaexon_t *metaexon = item->item;
	      num_left_breaks  += array_list_size(metaexon->left_breaks);
	      num_right_breaks += array_list_size(metaexon->right_breaks);
	      num_meta++;	      
	    }
	  }
	
	  //printf("Pack meta %i, %i, %i...\n", num_meta, num_left_breaks, num_right_breaks);
	  unsigned long *list_meta = (unsigned long *)malloc(sizeof(unsigned long)*num_meta*5);
	  unsigned long *list_left_breaks  = (unsigned long *)malloc(sizeof(unsigned long)*num_left_breaks*2);
	  unsigned long *list_right_breaks = (unsigned long *)malloc(sizeof(unsigned long)*num_right_breaks*2);
	  
	  MPI_metaexon_package(metaexons, num_meta, list_meta, 
			       num_left_breaks, list_left_breaks, 
			       num_right_breaks, list_right_breaks);

	  MPI_Send(&num_sj, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);	  
	  MPI_Send(&num_meta, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(&num_left_breaks, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(&num_right_breaks, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  //printf("[%i]SEND VALUES : %lu, %lu, %lu, %lu\n", rank, num_sj, num_meta, num_left_breaks, num_right_breaks);	  

	  MPI_Send(list_sj, num_sj*5, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(list_meta, num_meta*5, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(list_left_breaks, num_left_breaks*2, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(list_right_breaks, num_right_breaks*2, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);

	  free(list_sj);
	  free(list_meta);
	  free(list_left_breaks);
	  free(list_right_breaks);
	
	}
	
	if (rank == recv) {
	  //Recv AVLs
	  unsigned long recv_num_sj;	    
	  unsigned long recv_num_meta;
	  unsigned long recv_num_left_breaks;
	  unsigned long recv_num_right_breaks;
	
	  //printf("RECV: %i, %i, %i, %i, %i\n", recv_list_meta[recv_num_meta], recv_list_meta[recv_num_meta + 1], recv_list_meta[recv_num_meta + 2], recv_list_meta[recv_num_meta + 3], recv_list_meta[recv_num_meta + 4]);	  
	
	  MPI_Recv(&recv_num_sj, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);	
	  MPI_Recv(&recv_num_meta, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  MPI_Recv(&recv_num_left_breaks, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  MPI_Recv(&recv_num_right_breaks, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  //printf("[%i]RECV VALUES : %lu, %lu, %lu, %lu\n", rank, recv_num_sj, recv_num_meta, recv_num_left_breaks,recv_num_right_breaks);	  

	  unsigned long *recv_list_sj = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_sj*5);
	  unsigned long *recv_list_meta = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_meta*5);
	  unsigned long *recv_list_left_breaks = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_left_breaks*2);
	  unsigned long *recv_list_right_breaks = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_right_breaks*2);

	  MPI_Recv(recv_list_sj, recv_num_sj*5, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);	  	
	  MPI_Recv(recv_list_meta, recv_num_meta*5, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  MPI_Recv(recv_list_left_breaks, recv_num_left_breaks*2, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  MPI_Recv(recv_list_right_breaks, recv_num_right_breaks*2, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);


	  //Merge AVLs
	  avl_node_t *node_avl_start, *node_avl_end;
	  MPI_splice_t MPI_sj;
	  for (size_t sj = 0; sj < recv_num_sj; sj++) {
	    MPI_sj.start          = recv_list_sj[sj];
	    MPI_sj.end            = recv_list_sj[recv_num_sj + sj];
	    MPI_sj.reads_number   = recv_list_sj[recv_num_sj*2 + sj];
	    MPI_sj.strand         = recv_list_sj[recv_num_sj*3 + sj];
	    MPI_sj.chr            = recv_list_sj[recv_num_sj*4 + sj];
	  
	    //printf("\t[%i] %i/%i Merge (%i) %i: %lu - %lu\n", rank, sj, recv_num_sj, MPI_sj.strand, MPI_sj.chr, MPI_sj.start, MPI_sj.end);
	    allocate_start_node(MPI_sj.chr - 1,
				MPI_sj.strand,
				MPI_sj.start,
				MPI_sj.end,
				MPI_sj.start,
				MPI_sj.end,
				FROM_READ,
				1,
				"",
				&node_avl_start,
				&node_avl_end,
				avls_list);
	    //Refresh MPI_sj.reads_number in nodes
	    //end_data_t *end_data = node_avl_start->data;	    
	  }
	  //Merge Metaexon
	  //sleep(1);
	
	  //printf("[%i]Merge %i META:\n", rank, recv_num_meta);
	  size_t l_pos = 0, r_pos = 0;
	  MPI_metaexon_t mpi_meta;
	  for (size_t m = 0; m < recv_num_meta; m++) {
	    mpi_meta.chromosome = recv_list_meta[m];
	    mpi_meta.start      = recv_list_meta[recv_num_meta   + m];
	    mpi_meta.end        = recv_list_meta[recv_num_meta*2 + m];
	    mpi_meta.n_starts   = recv_list_meta[recv_num_meta*3 + m];
	    mpi_meta.n_ends     = recv_list_meta[recv_num_meta*4 + m];

	    if (!mpi_meta.n_starts && !mpi_meta.n_ends) {
	      //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends)
	      metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			      METAEXON_LEFT_END, NULL, metaexons);
	    } 
	    
	    //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
	    MPI_breaks_t bk;
	    for (size_t s = 0; s < mpi_meta.n_starts; s++) {
	      bk.pos    = recv_list_right_breaks[r_pos];
	      bk.strand = recv_list_right_breaks[recv_num_right_breaks + r_pos];
	      r_pos++;
	      //MPI_breaks_t bk = MPI_breaks[se_pos++];
	      //printf("\tM-STARTS-(%i)%lu\n", bk.strand, bk.pos);
	      cp_avltree *avl = avls_list->avls[bk.strand][mpi_meta.chromosome].avl;
	      avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	      assert(avl_node);
	      //if (!avl_node) {
	      //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
	      //printf("%i:%i:%i\n", m, r_pos - 1, bk.pos);
	      //exit(-1);
	      //}
	      metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			      METAEXON_RIGHT_END, avl_node, metaexons);
	    }
	    for (size_t s = 0; s < mpi_meta.n_ends; s++) {
	      bk.pos    = recv_list_left_breaks[l_pos];
	      bk.strand = recv_list_left_breaks[recv_num_left_breaks + l_pos];
	      l_pos++;
	      //MPI_breaks_t bk = MPI_breaks[se_pos++];
	      //printf("\t[%i]M-ENDS-(%i)%lu\n", rank, bk.strand, bk.pos);
	      cp_avltree *avl = avls_list->ends_avls[bk.strand][mpi_meta.chromosome].avl;
	      avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	      assert(avl_node);
	      metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			      METAEXON_LEFT_END, avl_node, metaexons);
	    }
	  }
	
	  free(recv_list_sj);
	  free(recv_list_meta);
	  free(recv_list_left_breaks);
	  free(recv_list_right_breaks);

	}
	send -= dec;
	recv -= dec;
	if (send <= 0) { send = -1; }
      }
      dec   *= 2;
      dec_s *= 2;
      dec_r *= 2;
    }
  
    //printf("[%i]MERGE NODO 0 WORKFLOW 1 END...BROADCAST TO OTHER NODES\n", rank);

    //int num_sj;
    unsigned long recv_num_sj, recv_n_metaexons, recv_n_starts_ends;
    //MPI_splice_t *MPI_splice;
    //MPI_metaexon_package_t *mpi_meta_pack;
    //MPI_metaexon_t *MPI_metaexon;
    //MPI_breaks_t *MPI_breaks;
  
    unsigned long *list_sj;
    unsigned long *list_meta;
    unsigned long *list_left_breaks;
    unsigned long *list_right_breaks;
    
    unsigned long num_meta = 0;
    unsigned long num_left_breaks = 0, num_right_breaks = 0;    
    
    if (rank == 0) {
      //PACK AVL
      num_sj = get_total_items(num_chromosomes, avls_list);
      list_sj = (unsigned long *)malloc(sizeof(unsigned long)*num_sj*5); 	  
      MPI_avl_package(num_chromosomes, avls_list, num_sj, list_sj);	  
      
      //PACK METAEXON
      for (int c = 0; c < metaexons->num_chromosomes; c++) {
	linked_list_t *metaexon_list = metaexons->metaexons_list[c];
	linked_list_item_t *item;
	for (item = metaexon_list->first; item != NULL; item = item->next) {
	  metaexon_t *metaexon = item->item;
	  num_left_breaks  += array_list_size(metaexon->left_breaks);
	  num_right_breaks += array_list_size(metaexon->right_breaks);
	  num_meta++;	      
	}
      }
      
      list_meta = (unsigned long *)malloc(sizeof(unsigned long)*num_meta*5);
      list_left_breaks  = (unsigned long *)malloc(sizeof(unsigned long)*num_left_breaks*2);
      list_right_breaks = (unsigned long *)malloc(sizeof(unsigned long)*num_right_breaks*2);

      MPI_metaexon_package(metaexons, num_meta, list_meta, 
			   num_left_breaks, list_left_breaks, 
			   num_right_breaks, list_right_breaks);

    }
    
    //printf("[%i]Bcast Send\n", rank);
    
    MPI_Bcast(&num_sj, 1, MPI_UNSIGNED_LONG, 0, new_comm);
    MPI_Bcast(&num_meta, 1, MPI_UNSIGNED_LONG, 0, new_comm);
    MPI_Bcast(&num_left_breaks, 1, MPI_UNSIGNED_LONG, 0, new_comm);
    MPI_Bcast(&num_right_breaks, 1, MPI_UNSIGNED_LONG, 0, new_comm);
    //printf("[%i]Pack meta %i %i, %i, %i...\n", rank, num_sj, num_meta, num_left_breaks, num_right_breaks);
    
    //size_t xn = 100000000;
    size_t xn = num_sj*5 + num_meta*5 + num_left_breaks*2 + num_right_breaks*2;
    unsigned long *pack_send = (unsigned long *)malloc(sizeof(unsigned long)*xn);

    if (rank == 0) {
      memcpy(pack_send, list_sj, sizeof(unsigned long)*num_sj*5);
      memcpy(&pack_send[num_sj*5], list_meta, sizeof(unsigned long)*num_meta*5);
      memcpy(&pack_send[num_sj*5 + num_meta*5], list_left_breaks, sizeof(unsigned long)*num_left_breaks*2);
      memcpy(&pack_send[num_sj*5 + num_meta*5 + num_left_breaks*2], list_right_breaks, sizeof(unsigned long)*num_right_breaks*2);

      free(list_sj);
      free(list_meta);
      free(list_left_breaks);
      free(list_right_breaks);
    }

    MPI_Bcast(pack_send, xn, MPI_UNSIGNED_LONG, 0, new_comm);

    //printf("[%i] RECV2: %lu, %lu, %lu \n", rank, pack_send[0], pack_send[1], pack_send[num_meta]);

    list_sj           = pack_send;
    list_meta         = &pack_send[num_sj*5];
    list_left_breaks  = &pack_send[num_sj*5 + num_meta*5];
    list_right_breaks = &pack_send[num_sj*5 + num_meta*5 + num_left_breaks*2];
  
    //  exit(-1);

    //printf("[%i]Merge 2 %i:\n", rank, num_sj);

    //Merge AVLs
    avl_node_t *node_avl_start, *node_avl_end;
    MPI_splice_t MPI_sj;

    for (size_t sj = 0; sj < num_sj; sj++) {
      MPI_sj.start          = list_sj[sj];
      MPI_sj.end            = list_sj[num_sj + sj];
      MPI_sj.reads_number   = list_sj[num_sj*2 + sj];
      MPI_sj.strand         = list_sj[num_sj*3 + sj];
      MPI_sj.chr            = list_sj[num_sj*4 + sj];

      //printf("[%i]Merge 2 (%i) %i: %lu - %lu\n", rank, MPI_sj.strand, MPI_sj.chr, MPI_sj.start, MPI_sj.end);
      allocate_start_node(MPI_sj.chr - 1,
			  MPI_sj.strand,
			  MPI_sj.start,
			  MPI_sj.end,
			  MPI_sj.start,
			  MPI_sj.end,
			  FROM_READ,
			  1,
			  "",
			  &node_avl_start,
			  &node_avl_end,
			  avls_list);
      //Refresh MPI_sj.reads_number in nodes
      end_data_t *end_data = node_avl_start->data;	    
    }
    
    //printf("[%i]Merge Meta %i\n", rank, num_meta);
    
    //Merge Metaexon
    size_t se_pos = 0;
    size_t r_pos = 0, l_pos = 0;
    MPI_metaexon_t mpi_meta;
    
    for (size_t m = 0; m < num_meta; m++) {
      mpi_meta.chromosome = list_meta[m];
      mpi_meta.start      = list_meta[num_meta   + m];
      mpi_meta.end        = list_meta[num_meta*2 + m];
      mpi_meta.n_starts   = list_meta[num_meta*3 + m];
      mpi_meta.n_ends     = list_meta[num_meta*4 + m];

      if (!mpi_meta.n_starts && !mpi_meta.n_ends) {
	//printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends)
	metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			METAEXON_LEFT_END, NULL, metaexons);
      } 
      
      //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
      MPI_breaks_t bk;
      for (int s = 0; s < mpi_meta.n_starts; s++) {
	bk.pos    = list_right_breaks[r_pos];
	bk.strand = list_right_breaks[num_right_breaks + r_pos];
	r_pos++;

	//printf("\tM-STARTS-(%i)%lu\n", bk.strand, bk.pos);
	cp_avltree *avl = avls_list->avls[bk.strand][mpi_meta.chromosome].avl;
	avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	//assert(avl_node);
	if (!avl_node) {
	  printf("ERROR: (%i|%i)%lu\n", m, s, bk.pos);
	  exit(-1);
	}
	metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			METAEXON_RIGHT_END, avl_node, metaexons);
      }
      for (int s = 0; s < mpi_meta.n_ends; s++) {
	bk.pos    = list_left_breaks[l_pos];
	bk.strand = list_left_breaks[num_left_breaks + l_pos];
	l_pos++;

	//printf("\t[%i]M-ENDS-(%i)%lu\n", rank, bk.strand, bk.pos);
	cp_avltree *avl = avls_list->ends_avls[bk.strand][mpi_meta.chromosome].avl;
	avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	assert(avl_node);
	metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			METAEXON_LEFT_END, avl_node, metaexons);
      }
    }
    
    //printf("Frees merge\n");
    free(pack_send);

    //====                E N D    M E R G E                ====//
    //==========================================================//

    t_merge_2 = 0;
    stop_timer(start_BWT_2, stop_BWT_2, t_merge_2);
    
    //if (rank == 0) {
    //printf("Merge 2 time: %0.2f\n", t_BWT_2 / 1000000);
    //}
    
    //End Writer
    char *end = (char *)malloc(sizeof(char)*MAX_BUFFER);
    strcpy(end, "END\0");    
    //printf("[%i]SEND  MPI END\n", rank);
    MPI_Send(end, MAX_BUFFER, MPI_CHAR, numprocs - 1, 1, MPI_COMM_WORLD);
    free(end);

    //printf("\t\t[%i]Task = %0.2fs, Send = %0.2fs\n", rank, t_task / 1000000, t_send / 1000000);
    
  } else {
    
    double t_write = 0;
    struct timeval start_task, stop_task;
    
    //printf("[%i]MPI THREAD WRITER START\n", rank);
    
    //=============================//
    //= W R I T E R - T H R E A D =//
    //=============================//
    //fd_out = fopen(output_filename, "a");
    
    char *recv_buff = (char *)malloc(sizeof(char)*MAX_BUFFER);
    int numprocs_end = 0;
    while (1) {
      MPI_Recv(recv_buff, MAX_BUFFER, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
      if (strcmp("END", recv_buff) == 0) {
	numprocs_end++;
	if (numprocs_end == numprocs - 1) {
	  break;
	}
      } else {
	fwrite(recv_buff, strlen(recv_buff), sizeof(char), fd_out);
      }
    }       
  }


  MPI_Barrier(MPI_COMM_WORLD);

  /*
  if (rank == 0) {
    stop_timer(wx_start, wx_stop, wx_time);
    printf("================================= WORKFLOW 2 TIME %0.2fs ========================\n", wx_time / 1000000);
    wx_time = 0;
    start_timer(wx_start);
  }
  */
 fast_mode_skip:
  
  reads_wf2 = reads_wf;
  reads_wf  = 0;
  wf_input_file_t *wf_input_file;
  workflow_t *wf_hc;
    
  workflow_SA_t *wf_last_SA;
  //== W3 ==//
  if (rank != numprocs - 1) {
    MPI_data.queue_tasks = queue_hc_tasks;
    if (options->fast_mode) {
      //Create and initialize second workflow
      wf_last_SA = workflow_SA_new();
      workflow_stage_function_SA_t stage_functions_last[] = {sa_rna_mapper_last};
      char *stage_labels_last[] = {"SA mapper last stage"};
      workflow_set_stages_SA(1, (workflow_stage_function_SA_t *)&stage_functions_last, stage_labels_last, wf_last_SA);
      workflow_set_producer_SA((workflow_producer_function_SA_t *)sa_alignments_reader_rna, "FastQ reader", wf_last_SA);
      workflow_set_consumer_SA((workflow_consumer_function_SA_t *)SA_MPI_consumer, "SAM writer", wf_last_SA);
    } else {
      wf_input_file = wf_input_file_new(NULL, batch);     
      wf_hc = workflow_new();
      workflow_stage_function_t stage_functions_hc[] = {rna_last_hc_stage, post_pair_stage};
      char *stage_labels_hc[] = {"RNA HARD CLIPPINGS", "POST PAIR"};
      workflow_set_stages(2, (workflow_stage_function_t *)&stage_functions_hc, stage_labels_hc, wf_hc);
      workflow_set_producer((workflow_producer_function_t *)file_reader_2, "Buffer reader", wf_hc);
      workflow_set_consumer((workflow_consumer_function_t *)MPI_output, "SAM writer", wf_hc);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank != numprocs - 1) {
    //=============================//
    //= W O R K E R - T H R E A D =//
    //=============================//
    //char *file;
    //while (1) {
    //if (MPI_data.queue_tasks->size) {
    //file = (char *)linked_list_remove_first(MPI_data.queue_tasks);
    //} else {
    //file = NULL;
    //}    
    //if (!file) { printf("END\n"); break; }      
    
    char buffer_hc_path[1024];
    sprintf(buffer_hc_path, "%s/buffer_hc.%i", TMP_PATH_BUFFERS, rank);      
    f_hc = fopen(buffer_hc_path, "r+b");
    if (f_hc == NULL) {
      printf("Error opening file '%s' \n", buffer_hc_path);
      exit(-1);
    }
      
    sw_input.f_hc = f_hc;

    if (options->fast_mode) {
      //sa_wf_input->file = f_hc;
      sa_rna.file1 = f_hc;
      struct timeval start_sa, stop_sa;
      double time_sa = 0;
      start_timer(start_sa);
      workflow_run_with_SA(options->num_cpu_threads, sa_wf_input, wf_last_SA);
      stop_timer(start_sa, stop_sa, time_sa);
      //printf("[%i]W3: Process SA file with %i threads in %0.2f(s)\n", rank, options->num_cpu_threads, time_sa / 1000000);	
      wf_last_SA->completed_producer = 0;
    } else {
      wf_input_file->file = f_hc;
      workflow_run_with(options->num_cpu_threads, wf_input_file, wf_hc);
      //workflow_display_timing(wf_hc);
      wf_hc->completed_producer = 0;
    }
    
    fclose(f_hc);
    
    batch_out_t *b = array_list_get(1, data_out->buffer);    
    MPI_Send(b->data, MAX_BUFFER, MPI_CHAR, numprocs - 1, 1, MPI_COMM_WORLD);
    b->len = 0;
    
    //End Writer
    char *end = (char *)malloc(sizeof(char)*MAX_BUFFER);
    strcpy(end, "END\0");
    MPI_Send(end, MAX_BUFFER, MPI_CHAR, numprocs - 1, 1, MPI_COMM_WORLD);
    
  } else {

    double t_write = 0;
    struct timeval start_task, stop_task;
    
    //printf("[%i]MPI THREAD WRITER START\n", rank);
    
    //=============================//
    //= W R I T E R - T H R E A D =//
    //=============================//
    //fd_out = fopen(output_filename, "a");
    
    char *recv_buff = (char *)malloc(sizeof(char)*MAX_BUFFER);
    int numprocs_end = 0;
    while (1) {
      MPI_Recv(recv_buff, MAX_BUFFER, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
      if (strcmp("END", recv_buff) == 0) {
	numprocs_end++;
	if (numprocs_end == numprocs - 1) {
	  break;
	}
      } else {
	fwrite(recv_buff, strlen(recv_buff), sizeof(char), fd_out);
      }    
    }
    
  }

  MPI_Barrier(MPI_COMM_WORLD);
  /*  
  if (rank == 0) {
    wx_time = 0;
    stop_timer(wx_start, wx_stop, wx_time);
    printf("WORKFLOW 3 TIME %0.2fs\n", wx_time / 1000000);
    start_timer(wx_start);
  }
  */
  reads_wf3 = reads_wf;
  reads_wf  = 0;
   
  //============================================================//
  //= M E R G E    M E T A E X O N - A V L    T O    W R I T E =//
  //============================================================//
  //== Package structures ==//


  //////////////////////////////////////////////////////////////////
  MPI_Barrier(MPI_COMM_WORLD);
  /*
  if (rank != numprocs - 1) {
    unsigned long num_sj = get_total_items(num_chromosomes, avls_list);
    //SEND METAEXON
    unsigned long num_meta = 0;
    unsigned long num_left_breaks = 0, num_right_breaks = 0;
    for (int c = 0; c < metaexons->num_chromosomes; c++) {
      linked_list_t *metaexon_list = metaexons->metaexons_list[c];
      linked_list_item_t *item;
      for (item = metaexon_list->first; item != NULL; item = item->next) {
	metaexon_t *metaexon = item->item;
	num_left_breaks  += array_list_size(metaexon->left_breaks);
	num_right_breaks += array_list_size(metaexon->right_breaks);
	num_meta++;	      
      }
    }
    printf("3.[%i]num_sj = %lu, num_meta = %lu\n", rank, num_sj, num_meta);
    msj_3 = num_sj;
    mmeta_3 = num_meta;
  }
*/
  //////////////////////////////////////////////////////////////////

  struct timeval start_BWT_3, stop_BWT_3;
  double t_merge_3;

  start_timer(start_BWT_3);  

  if (rank != numprocs - 1) {

    unsigned long num_sj;        
    int numprocs_merge = numprocs - 1;
    int dec_s = 1;
    int dec_r = 2;
    int send, recv;
    int dec = 2;
    
    while (1) {
      send = numprocs_merge - dec_s;
      recv = numprocs_merge - dec_r;
      if (send > 0 && recv < 0)
	recv = 0;
      if (send <= 0 && recv <= 0)
	break;
      
      for (int r = numprocs_merge - 1; r >= 0; r--) {
	if (rank == send) {
	  //SEND AVL
	  num_sj = get_total_items(num_chromosomes, avls_list);
	  unsigned long  *list_sj = (unsigned long *)malloc(sizeof(unsigned long)*num_sj*5); 	  
	  MPI_avl_package(num_chromosomes, avls_list, num_sj, list_sj);	  

	  //SEND METAEXON
	  unsigned long num_meta = 0;
	  unsigned long num_left_breaks = 0, num_right_breaks = 0;
	  for (int c = 0; c < metaexons->num_chromosomes; c++) {
	    linked_list_t *metaexon_list = metaexons->metaexons_list[c];
	    linked_list_item_t *item;
	    for (item = metaexon_list->first; item != NULL; item = item->next) {
	      metaexon_t *metaexon = item->item;
	      num_left_breaks  += array_list_size(metaexon->left_breaks);
	      num_right_breaks += array_list_size(metaexon->right_breaks);
	      num_meta++;	      
	    }
	  }
	
	  //printf("Pack meta %i, %i, %i...\n", num_meta, num_left_breaks, num_right_breaks);
	  unsigned long *list_meta = (unsigned long *)malloc(sizeof(unsigned long)*num_meta*5);
	  unsigned long *list_left_breaks  = (unsigned long *)malloc(sizeof(unsigned long)*num_left_breaks*2);
	  unsigned long *list_right_breaks = (unsigned long *)malloc(sizeof(unsigned long)*num_right_breaks*2);
	  
	  MPI_metaexon_package(metaexons, num_meta, list_meta, 
			       num_left_breaks, list_left_breaks, 
			       num_right_breaks, list_right_breaks);

	  //printf("RECV-RANK0: %i, %i, %i, %i, %i\n", list_meta[num_meta], list_meta[num_meta + 1], list_meta[num_meta + 2], list_meta[num_meta + 3], list_meta[num_meta + 4]);

	  MPI_Send(&num_sj, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);	  
	  MPI_Send(&num_meta, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(&num_left_breaks, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(&num_right_breaks, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  //printf("[%i]SEND VALUES : %lu, %lu, %lu, %lu\n", rank, num_sj, num_meta, num_left_breaks, num_right_breaks);	  

	  MPI_Send(list_sj, num_sj*5, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(list_meta, num_meta*5, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(list_left_breaks, num_left_breaks*2, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  MPI_Send(list_right_breaks, num_right_breaks*2, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);

	  //printf("[%i]SEND LISTS : %lu, %lu, %lu, %lu\n", rank, list_sj[num_sj], list_meta[num_meta], list_sj[num_sj + 1], list_meta[num_meta + 1]);
	  //printf("SEND: %i, %i, %i, %i, %i\n", list_right_breaks[0], list_right_breaks[1], list_right_breaks[2], list_right_breaks[3], list_right_breaks[4]);
	
	  free(list_sj);
	  free(list_meta);
	  free(list_left_breaks);
	  free(list_right_breaks);
	
	}
	
	if (rank == recv) {
	  //Recv AVLs
	  unsigned long recv_num_sj;	    
	  unsigned long recv_num_meta;
	  unsigned long recv_num_left_breaks;
	  unsigned long recv_num_right_breaks;
	
	  //printf("RECV: %i, %i, %i, %i, %i\n", recv_list_meta[recv_num_meta], recv_list_meta[recv_num_meta + 1], recv_list_meta[recv_num_meta + 2], recv_list_meta[recv_num_meta + 3], recv_list_meta[recv_num_meta + 4]);	  
	
	  MPI_Recv(&recv_num_sj, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);	
	  MPI_Recv(&recv_num_meta, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  MPI_Recv(&recv_num_left_breaks, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  MPI_Recv(&recv_num_right_breaks, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  //printf("[%i]RECV VALUES : %lu, %lu, %lu, %lu\n", rank, recv_num_sj, recv_num_meta, recv_num_left_breaks,recv_num_right_breaks);	  

	  unsigned long *recv_list_sj = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_sj*5);
	  unsigned long *recv_list_meta = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_meta*5);
	  unsigned long *recv_list_left_breaks = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_left_breaks*2);
	  unsigned long *recv_list_right_breaks = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_right_breaks*2);

	  MPI_Recv(recv_list_sj, recv_num_sj*5, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);	  	
	  MPI_Recv(recv_list_meta, recv_num_meta*5, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  MPI_Recv(recv_list_left_breaks, recv_num_left_breaks*2, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  MPI_Recv(recv_list_right_breaks, recv_num_right_breaks*2, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);


	  //Merge AVLs
	  avl_node_t *node_avl_start, *node_avl_end;
	  MPI_splice_t MPI_sj;
	  for (size_t sj = 0; sj < recv_num_sj; sj++) {
	    MPI_sj.start          = recv_list_sj[sj];
	    MPI_sj.end            = recv_list_sj[recv_num_sj + sj];
	    MPI_sj.reads_number   = recv_list_sj[recv_num_sj*2 + sj];
	    MPI_sj.strand         = recv_list_sj[recv_num_sj*3 + sj];
	    MPI_sj.chr            = recv_list_sj[recv_num_sj*4 + sj];

	    
	    //printf("\t[%i] %i/%i Merge (%i) %i: %lu - %lu\n", rank, sj, recv_num_sj, MPI_sj.strand, MPI_sj.chr, MPI_sj.start, MPI_sj.end);
	    allocate_start_node(MPI_sj.chr - 1,
				MPI_sj.strand,
				MPI_sj.start,
				MPI_sj.end,
				MPI_sj.start,
				MPI_sj.end,
				FROM_READ,
				1,
				"",
				&node_avl_start,
				&node_avl_end,
				avls_list);
	    //Refresh MPI_sj.reads_number in nodes
	    //end_data_t *end_data = node_avl_start->data;	    
	  }
	  //Merge Metaexon
	  //sleep(1);
	
	  //printf("[%i]Merge %i META:\n", rank, recv_num_meta);
	  size_t l_pos = 0, r_pos = 0;
	  MPI_metaexon_t mpi_meta;
	  for (size_t m = 0; m < recv_num_meta; m++) {
	    mpi_meta.chromosome = recv_list_meta[m];
	    mpi_meta.start      = recv_list_meta[recv_num_meta   + m];
	    mpi_meta.end        = recv_list_meta[recv_num_meta*2 + m];
	    mpi_meta.n_starts   = recv_list_meta[recv_num_meta*3 + m];
	    mpi_meta.n_ends     = recv_list_meta[recv_num_meta*4 + m];

	    if (!mpi_meta.n_starts && !mpi_meta.n_ends) {
	      //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends)
	      metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			      METAEXON_LEFT_END, NULL, metaexons);
	    } 

	    //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
	    MPI_breaks_t bk;
	    for (size_t s = 0; s < mpi_meta.n_starts; s++) {
	      bk.pos    = recv_list_right_breaks[r_pos];
	      bk.strand = recv_list_right_breaks[recv_num_right_breaks + r_pos];
	      r_pos++;
	      //MPI_breaks_t bk = MPI_breaks[se_pos++];
	      //printf("\tM-STARTS-(%i)%lu\n", bk.strand, bk.pos);
	      cp_avltree *avl = avls_list->avls[bk.strand][mpi_meta.chromosome].avl;
	      avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	      assert(avl_node);
	      //if (!avl_node) {
	      //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
	      //printf("%i:%i:%i\n", m, r_pos - 1, bk.pos);
	      //exit(-1);
	      //}
	      metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			      METAEXON_RIGHT_END, avl_node, metaexons);
	    }
	    for (size_t s = 0; s < mpi_meta.n_ends; s++) {
	      bk.pos    = recv_list_left_breaks[l_pos];
	      bk.strand = recv_list_left_breaks[recv_num_left_breaks + l_pos];
	      l_pos++;
	      //MPI_breaks_t bk = MPI_breaks[se_pos++];
	      //printf("\t[%i]M-ENDS-(%i)%lu\n", rank, bk.strand, bk.pos);
	      cp_avltree *avl = avls_list->ends_avls[bk.strand][mpi_meta.chromosome].avl;
	      avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	      assert(avl_node);
	      metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			      METAEXON_LEFT_END, avl_node, metaexons);
	    }
	  }
	
	  free(recv_list_sj);
	  free(recv_list_meta);
	  free(recv_list_left_breaks);
	  free(recv_list_right_breaks);

	}
	send -= dec;
	recv -= dec;
	if (send <= 0) { send = -1; }
      }
      dec   *= 2;
      dec_s *= 2;
      dec_r *= 2;

    }

    MPI_Barrier(new_comm);
    
    t_merge_3 = 0;
    stop_timer(start_BWT_3, stop_BWT_3, t_merge_3);
    
    //if (rank == 0) {
    //printf("Merge 3 time: %0.2f\n", t_BWT_3 / 1000000);
    //}

    //printf("[%i]MERGE NODO 0 WORKFLOW 1 END...BROADCAST TO OTHER NODES\n", rank);

    if (rank == 0) {
      write_chromosome_avls(NULL, exact_filename, num_chromosomes, avls_list);      
    }
    
    //
    // end of workflow management
    //--------------------------------------------------------------------------------------
    
    // free memory
    //workflow_free(wf);
    //wf_input_free(wf_input);
    //Write chromosome avls
    //batch_free(batch);

    //free(output_filename);
    //free(exact_filename);
  } else {
    fclose(fd_out);
  }


  MPI_Barrier(MPI_COMM_WORLD);
  /*
  if (rank == 0) {
    wx_time = 0;
    stop_timer(wx_start, wx_stop, wx_time);
    printf("WRITE SJ TIME %0.2fs\n", wx_time / 1000000);
  }
  */
  
  if (rank == 0) {
    //stop_timer(start_w, stop_w, t_m3);
    stop_timer(start_total, stop_total, t_total);
    //printf("\tSTOP MERGE -3- (%0.2fs)\n", t_m3 / 1000000);
  }

  if (rank != numprocs - 1) {
    // free memory
    if (bwt_index) { bwt_index_free(bwt_index); }
    metaexons_free(metaexons);
    
    if (genome)
      genome_free(genome);

  }
  /*
  if (rank == 0) {
    printf("MERGEs: %0.2f | %0.2f | %0.2f\n", rank, t_merge_1 / 1000000, t_merge_2 / 1000000, t_merge_3 / 1000000);
  }
  
  if (rank != numprocs -1) {
    printf("[%i] SJ:     %lu | %lu | %lu\n", rank, msj_1, msj_2, msj_3);
    printf("[%i] MERGE:  %lu | %lu | %lu\n", rank, mmeta_1, mmeta_2, mmeta_3);
  }
  */
  if (rank == 0) {
    printf("==================================================================\n");
    printf(" TOTAL TIME          : %0.2f(s)\n", t_total / 1000000);
    printf("==================================================================\n\n\n");
  }

  //printf("[%i]END\n", rank);

 exit:
  MPI_Finalize();
  
}




void rna_aligner_mpi_work_stealing(options_t *options, int argc, char *argv[]) {
  //MPI Process Files  in the nodes
  FILE* fd_out;
  int tag = 1;
  MPI_Status status;
  int numprocs;
  int rank;
  int len;  
  n_tasks_t ntasks;
  int th_steal = 0;
  int worker   = 1;
  int total_files;
  MPI_data_t MPI_data;
  
  MPI_data.block = 1;
  MPI_data.run = 1;
  pthread_mutex_init(&(MPI_data.queue_mutex), NULL);
  //Queue of tasks
  // ([0] Task for process this node) <--- [0]-[1]-[2]-[3] ---> ([3] Task for steal from thief node)
  MPI_data.queue_tasks = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);

  linked_list_t *queue_sa_tasks = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  linked_list_t *queue_hc_tasks = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  
  
  int total_tasks = 0;
  int init_id;
  
  argc += 1;
  argv -= 1;
  
  MPI_Init(&argc,&argv);
  int provided;
  //MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  //if (provided != MPI_THREAD_MULTIPLE) {
  //printf("Sorry, this MPI implementation does not support multiple threads\n");
  //MPI_Abort(MPI_COMM_WORLD, 1);
  //}
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  //================= S Y N C ==================//
  //unsigned char sync_data[numprocs];
  //memset(sync_data, 0, numprocs);
  //--------------------------------------------//

  
  //============== T I M I N G =================//
  double t_genomes, t_split, t_process, t_merge;               //t_file, t_file_write, t_total_file, t_loop;
  struct timeval start, stop, start_split, stop_split;         //start_file, stop_file, start_loop, stop_loop;

  double t_total = 0, t_w1 = 0, t_w2 = 0, t_w3 = 0, t_m1 = 0, t_m2 = 0, t_m3 = 0;

  struct timeval start_total, stop_total, start_w, stop_w;
  
  t_genomes = 0;
  t_split = 0;
  t_process = 0;
  t_merge = 0;

  double t_file, time_sp;
  struct timeval start_file, stop_file, start_sp, stop_sp, end_sp;

  t_file = 0;
  //-------------------------------------------//

  
  //========================= F I L E S    P A T H S ==========================//
  //End fill 
  int path_length = strlen(options->output_name);
  int prefix_length = 0;
  
  if (options->prefix_name) {
    prefix_length = strlen(options->prefix_name);
  }
  
  char *reads_results = (char *)calloc((60 + prefix_length), sizeof(char));
  char *exact_junctions = (char *)calloc((60 + prefix_length), sizeof(char));
  char *output_filename = (char *)calloc((path_length + prefix_length + 60), sizeof(char));
  char *exact_filename = (char *)calloc((path_length + prefix_length + 60), sizeof(char));
  //struct timeval time_genome_s, time_genome_e;
  //double time_genome = 0.0;
  metaexons_t *metaexons;
  genome_t *genome;  
  bwt_index_t *bwt_index = NULL;
  int num_chromosomes;
  // load index for dna/rna or for bisulfite case

  char filename_tab[strlen(options->bwt_dirname) + 1024];  
  sprintf(filename_tab, "%s/params.info", options->bwt_dirname);
  FILE *fd = fopen(filename_tab, "r");

  if (fd) { 
    options->fast_mode = 1; 
  } else { 
    options->fast_mode = 0; 
  }
  fclose(fd);
  
  if (options->fast_mode && !options->set_cal) {
    options->min_cal_size = 30;
  }

  if (options->fast_mode && options->set_cal) {
    if (options->min_cal_size < 30) {
      options->min_cal_size = 30;
    } else if(options->min_cal_size > 100)  {
      options->min_cal_size = 100;
    }
  }

  options->bam_format = 0;

  if (options->prefix_name) {
    strcat(reads_results, "/");
    strcat(reads_results, options->prefix_name);
    strcat(reads_results, "_");
    strcat(reads_results, OUTPUT_FILENAME);
    strcat(reads_results, ".sam");  
    strcat(exact_junctions, "/");
    strcat(exact_junctions, options->prefix_name);
    strcat(exact_junctions, "_exact_junctions.bed");
  } else {
    strcat(reads_results, "/");    
    strcat(reads_results, OUTPUT_FILENAME);
    strcat(reads_results, ".sam");      
    strcat(exact_junctions, "/exact_junctions.bed");
  } 

  strcat(output_filename, options->output_name);
  strcat(output_filename, reads_results);
  free(reads_results);

  strcat(exact_filename, options->output_name);
  strcat(exact_filename, exact_junctions);
  free(exact_junctions);

  char *file1 = options->in_filename;

  //printf("START: MPI-Worker-Thread Rank-%i\n", rank);
  //printf("MPI Hpg-aligner START in rank = %i, numprocs = %i!\n", rank, numprocs);

  //const int READS_SPLIT = 100000;
  const int READS_SPLIT = 100000;
  //const int READS_SPLIT = 2;
  
  if (rank == 0) {
    //printf("START MPI WITH %i PROCESS!\n", numprocs);
    start_timer(start);
    //==== S P L I T    I N P U T    F I L E ====//
    //TODO: Implement
    //create_directory(TMP_PATH_OUTPUTS);
    create_directory(TMP_PATH_BUFFERS);

    //#pragma omp parallel sections num_threads(2)
    //{
    //#pragma omp section
    //{
    //struct timeval start_sp, end_sp;
    //double time_sp;
    start_timer(start_split);
    //===================================================================//
    //                               O M P                               //
    //===================================================================//
    //                S P L I T    I N P U T    F I L E                  //
    //===================================================================//
    int total_files = split_input_file(file1, TMP_PATH_FILES, READS_SPLIT);
    //===================================================================//
    //                                                                   //
    //===================================================================//
    stop_timer(start_split, stop_split, t_split);
    printf("\t[%i]Split Input File Done! (%0.2fs)\n", rank, t_split / 1000000);
    
    //==== D E A L    W O R K  ====//
    int tasks_per_node = total_files / (numprocs - 1);
    int tasks_inc = total_files % (numprocs - 1);
    int tasks_index = 0;
    int actual_id;	
    init_id = 0;
    total_tasks = tasks_per_node;
    if (tasks_inc > 0) { total_tasks++; }    
    actual_id = total_tasks;	
    for (int se = 1; se < numprocs - 1; se++) {
      int send_tasks = tasks_per_node;
      if (se <= tasks_inc - 1) { send_tasks++; }      
      ntasks.num_tasks = send_tasks;
      ntasks.first_id  = actual_id;
      MPI_Send(&ntasks, sizeof(n_tasks_t), MPI_CHAR, se, 1, MPI_COMM_WORLD);      
      actual_id += send_tasks;
    }
    
    //}
    //#pragma omp section
	  //{
    //struct timeval start_sp, end_sp;
    //double time_sp;
    start_timer(start_sp);
    //===================================================================//
    //                L O A D    S T R U C T U R E S                     //
    //===================================================================//
    //////////////// LOAD BWT INDEX //////////////////////    
    bwt_index = bwt_index_new(options->bwt_dirname, false);
    ////////////////     GENOME     //////////////////////    
    genome = genome_new("dna_compression.bin", options->bwt_dirname, BWT_MODE);  
    num_chromosomes = genome->num_chromosomes;
    // Metaexons structure
    metaexons = metaexons_new(genome->num_chromosomes, 
			      genome->chr_size);
    //===================================================================//
    //                                                                   //
    //===================================================================//
    stop_timer(start_sp, end_sp, time_sp);
    printf("\t[%i]Genome Load Done! (%0.2fs)\n", rank, time_sp / 1000000);
	
    FILE* fd_p = fopen(output_filename, "w");
    write_sam_header_BWT(genome, (FILE *) fd_p);
    fclose(fd_p);
    //}
    //}
  } else if (rank != numprocs - 1) {
    //struct timeval start_sp, end_sp;
    //double time_sp;
    start_timer(start_sp);
    //===================================================================//
    //                L O A D    S T R U C T U R E S                     //
    //===================================================================//
    //////////////// LOAD BWT INDEX //////////////////////    
    bwt_index = bwt_index_new(options->bwt_dirname, false);
    ////////////////     GENOME     //////////////////////    
    genome = genome_new("dna_compression.bin", options->bwt_dirname, BWT_MODE);  
    num_chromosomes = genome->num_chromosomes;
    // Metaexons structure
    metaexons = metaexons_new(genome->num_chromosomes, 
			      genome->chr_size);
    //===================================================================//
    //                                                                   //
    //===================================================================//
    stop_timer(start_sp, end_sp, time_sp);
    printf("\t[%i]Genome Load Done! (%0.2fs)\n", rank, time_sp / 1000000);
    
    MPI_Recv(&ntasks, sizeof(n_tasks_t), MPI_CHAR, 0, 1, MPI_COMM_WORLD, &status);
    init_id   = ntasks.first_id;
    total_tasks = ntasks.num_tasks;    
  }

  bwt_optarg_t *bwt_optarg;
  cal_optarg_t *cal_optarg;
  pair_mng_t *pair_mng;
  report_optarg_t *report_optarg;
  fastq_batch_reader_input_t reader_input;
  linked_list_t *alignments_list;
  sw_optarg_t sw_optarg;  
  bwt_server_input_t bwt_input;
  region_seeker_input_t region_input;
  cal_seeker_input_t cal_input;
  pair_server_input_t pair_input;
  batch_writer_input_t writer_input;
  sw_server_input_t sw_input;
  avls_list_t* avls_list;
  
  workflow_t *wf;    
  wf_input_t *wf_input;

  FILE *f_sa, *f_hc;
  batch_t *batch;
  batch_buffer_t *data_out;
  
  if (rank != numprocs - 1) {
    //==== T A S K S    I N S E R T S ====//
    int item = init_id;
    for (int it = 0; it < total_tasks; it++) {
      char file_name[1024];
      sprintf(file_name, "%i.tmp", item);
      linked_list_insert_first(strdup(file_name), MPI_data.queue_tasks);
      linked_list_insert_first(strdup(file_name), queue_sa_tasks);
      linked_list_insert_first(strdup(file_name), queue_hc_tasks);
      item++;
    }

    printf("[%i]MODE: WORK-STEALING  - TASKS: %lu\n", rank, total_tasks);
	
    //=================== I N P U T    I N I T I A L I Z A T I O N S =========================//
    //                                                                                        //
    //========================================================================================// 
    //BWT parameters
    extern int min_intron, max_intron;
    min_intron = options->min_intron_length;
    max_intron = options->max_intron_length;
  
    bwt_optarg = bwt_optarg_new(1, 0,
				options->filter_read_mappings, 
				options->filter_seed_mappings);  
    // CAL parameters
    //printf("%i\n", options->min_cal_size);
    cal_optarg = cal_optarg_new(options->min_cal_size, 
				options->seeds_max_distance, 
				options->num_seeds, 
				options->min_num_seeds_in_cal,
				options->seed_size, 
				options->min_seed_size, 
				options->cal_seeker_errors, 
				options->max_intron_length, 
				options->min_intron_length);
    // paired mode parameters
    pair_mng = pair_mng_new(options->pair_mode, options->pair_min_distance, 
			    options->pair_max_distance, options->report_only_paired);
    // report parameters
    report_optarg = report_optarg_new(options->report_all,
				      options->report_n_best,
				      options->report_n_hits, 
				      options->report_only_paired,
				      options->report_best);

    //if (options->transcriptome_filename != NULL) {
    //printf("\nLoading transcriptome...\n");
    //load_transcriptome(options->transcriptome_filename, genome, avls_list, metaexons);
    //printf("Load done!\n");
    //}

    linked_list_t *buffer    = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);
    linked_list_t *buffer_hc = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);
    
    fastq_batch_reader_input_init(options->in_filename, options->in_filename2, 
				  options->pair_mode, options->batch_size, 
				  NULL, options->gzip, &reader_input);  


    alignments_list = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);
    linked_list_set_flag(options->pair_mode, alignments_list);

    sw_optarg_init(options->gap_open, options->gap_extend, 
		   options->match, options->mismatch, &sw_optarg);

    bwt_server_input_init(NULL,  0,  bwt_optarg, 
			  bwt_index, NULL,  0, 
			  NULL, metaexons, &sw_optarg, 
			  genome, &bwt_input);
      
    region_seeker_input_init(NULL, cal_optarg, 
			     bwt_optarg, bwt_index, NULL, 
			     0, 0, options->min_seed_padding_left, 
			     options->min_seed_padding_right, genome, metaexons, &region_input);
  
    cal_seeker_input_init(NULL, cal_optarg, NULL, 0, NULL, NULL, genome, 
			  bwt_optarg, bwt_index, metaexons, &cal_input);
  
    int pair_mode = pair_mng->pair_mode;

    avls_list = avls_list_new(num_chromosomes);
    sw_server_input_init(NULL, NULL, 0,  options->match,  
			 options->mismatch,  options->gap_open, options->gap_extend,  
			 options->min_score,  options->flank_length, genome,  
			 options->max_intron_length, options->min_intron_length,  
			 options->seeds_max_distance,  bwt_optarg, avls_list, 
			 cal_optarg, bwt_index, metaexons, buffer, buffer_hc, 
			 NULL, NULL, pair_mode, &sw_input);
    
    pair_server_input_init(pair_mng, report_optarg, NULL, NULL, NULL, &pair_input);

    batch_writer_input_init(output_filename,
			    exact_filename,
			    NULL,
			    //extend_filename, 
			    alignments_list, 
			    genome, 
			    &writer_input);

    writer_input.bam_format = options->bam_format;


    //clear_batch_buffer(batch_buffer_t *bb)

    data_out = new_batch_buffer();
  
    //init_data_output(&data_out, READS_SPLIT);
    batch = batch_new(&bwt_input, &region_input, &cal_input, 
		      &pair_input, NULL, &sw_input, &writer_input, RNA_MODE, NULL, data_out);

    //Parse input for more than one file
    char *file2 = options->in_filename2;

    fastq_batch_reader_input_init(file1, file2, 
				  options->pair_mode, options->batch_size, 
				  NULL, options->gzip, &reader_input);  

    fflush(stdout);
    
    wf_input = wf_input_new(&reader_input, batch);
  
    //create and initialize workflow
    wf = workflow_new();
    workflow_stage_function_t stage_functions[] = {bwt_stage, cal_stage, 
						   sw_stage, post_pair_stage, MPI_output};  
    char *stage_labels[] = {"BWT", "CAL", "SW", "POST PAIR", "MPI OUTPUT"};
    workflow_set_stages(5, (workflow_stage_function_t *)&stage_functions, stage_labels, wf);
    // optional producer and consumer functions
    workflow_set_producer((workflow_producer_function_t *)fastq_reader, "FastQ reader", wf);
  
    //if (options->bam_format) {
    //workflow_set_consumer((workflow_consumer_function_t *)bam_writer, "BAM writer", wf);
    //} else {
    workflow_set_consumer((workflow_consumer_function_t *)MPI_consumer, "SAM writer", wf);
    //workflow_set_consumer((workflow_consumer_function_t *)sam_writer, "SAM writer", wf);
    //}

    //========================================================================================//
    //                                                                                        //
    //========================================================================================// 
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0) {
    start_timer(start_total);
  }

  if (rank == 0) {
    printf("======================== WORKFLOW 1 ============================\n");
    start_timer(start_w);
  }

  double tot_t_file = 0;
  double tot_t_send = 0;

  double t_task = 0, t_send = 0;
  struct timeval start_task, stop_task, start_send, stop_send;
  int w1_tasks = 0;
  
  if (rank != numprocs - 1) {    
    size_t num_cpus = get_optimal_cpu_num_threads();
    //num_cpus -= 2;    
    if (rank == 0) {
      stop_timer(start, stop, t_genomes);
      start_timer(start);
    }
    
    //==== C R E A T E    P O L L I N G    T H R E A D ====//
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
    cpu_set_t cpu_set;
    CPU_ZERO(&cpu_set);
    CPU_SET(num_cpus - 1, &cpu_set);
    sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set), &cpu_set);
    
    int ret;
    pthread_t thread;
    if ((ret = pthread_create(&thread, &attr, mpi_thread_polling, (void *)&MPI_data))) {
      printf("ERROR; return code from pthread_create() is %d\n", ret);
      exit(-1);
    }
    
    //=============================//
    //= W O R K E R - T H R E A D =//
    //=============================//
    char *file;
    int victim;
    work_package_t work;
    
    //start_timer(start_loop);
    while (1) {
      //Process all tasks in the queue
      while (1) {
	pthread_mutex_lock(&(MPI_data.queue_mutex));
	if (MPI_data.queue_tasks->size) {
	  file = (char *)linked_list_remove_first(MPI_data.queue_tasks);
	} else {
	  file = NULL;
	}
	pthread_mutex_unlock(&(MPI_data.queue_mutex));
	if (!file) { break; }      

	start_timer(start_file);
	//====       W O R K E R    P R O C E S S    F I L E       ====//
	//=============================================================//
	char file_path[1024];
	char buffer_sa_path[1024];
	char buffer_hc_path[1024];
      
	sprintf(file_path, "%s/%s", TMP_PATH_FILES, file);
	sprintf(buffer_sa_path, "%s/buffer_sa.%s", TMP_PATH_BUFFERS, file);
	sprintf(buffer_hc_path, "%s/buffer_hc.%s", TMP_PATH_BUFFERS, file);
      
	//printf("rank.%i process file -> %s\n", rank, file_path);
	reader_input.fq_file1 = fastq_fopen(file_path);
      
	//Run workflow
	f_sa = fopen(buffer_sa_path, "w+b");
	if (f_sa == NULL) {
	  LOG_FATAL("Error opening file 'buffer_sa.tmp' \n");
	}
	sw_input.f_sa = f_sa;
	
	f_hc = fopen(buffer_hc_path, "w+b");
	if (f_hc == NULL) {
	  LOG_FATAL("Error opening file 'buffer_hc.tmp' \n");
	}

	start_timer(start_task);
	sw_input.f_hc = f_hc;
	t_file = 0;
	workflow_run_with(options->num_cpu_threads, wf_input, wf);
	w1_tasks++;
	tot_t_file += t_file;
	stop_timer(start_task, stop_task, t_task);
	
	//Send Batches to Write
	MPI_Request req;
	size_t len_send = (array_list_size(data_out->buffer) * MAX_BUFFER) / 1024; //KB
	len_send = len_send / 1024; //MB

	for (int i = 0; i < array_list_size(data_out->buffer); i++) {
	  batch_out_t *bo = array_list_get(i, data_out->buffer);
	  MPI_Send(bo->data, MAX_BUFFER, MPI_CHAR, numprocs - 1, 1, MPI_COMM_WORLD);
	  free_batch_out(bo);
	}
	array_list_clear(data_out->buffer, (void *)NULL);
	
	wf->completed_producer = 0;
	fastq_fclose(reader_input.fq_file1);
	fclose(f_sa);
	fclose(f_hc);
	
	//==== W O R K E R    E N D    P R O C E S S    F I L E ====//
	//==========================================================//
	stop_timer(start_file, stop_file, t_file);
	t_file = 0;
      }
      //Now I am idle. So, I send message to thread MPI comunicator
      work.source = rank;
      
      if (MPI_data.run) {
	work.mode = SAME_NODE;
	MPI_Send(&work, sizeof(work_package_t), MPI_CHAR, work.source, 1, MPI_COMM_WORLD);
      } else {
	work.mode = KILL_TH;
	MPI_Send(&work, sizeof(work_package_t), MPI_CHAR, work.source, 1, MPI_COMM_WORLD);
	break;
      }
      
      while (MPI_data.block) {}
      MPI_data.block = 1;

      //printf("\t[%i]Break block!\n", rank);
    }
    
    // free attribute and wait for the other threads
    void *pth_status;
    pthread_attr_destroy(&attr);
    
    if ((ret = pthread_join(thread, &pth_status))) {
      printf("ERROR; return code from pthread_join() is %d\n", ret);
      exit(-1);
    }
    
    //End Writer
    char *end = (char *)malloc(sizeof(char)*MAX_BUFFER);
    strcpy(end, "END\0");
    MPI_Send(end, MAX_BUFFER, MPI_CHAR, numprocs - 1, 1, MPI_COMM_WORLD);

    printf("[%i] Total Tasks %i process in %0.2fs (Time x task = %0.2fs), Total reads process %lu\n", rank, w1_tasks, t_task / 1000000, (t_task / 1000000) / w1_tasks, mpi_tot_reads);
    //printf("\t[%i] Total Tasks %i process in %0.2fs\n", rank, w1_tasks, t_task / 1000000);    

  } else {
    //=============================//
    //= W R I T E R - T H R E A D =//
    //=============================//
    fd_out = fopen(output_filename, "a");
    char *recv_buff = (char *)malloc(sizeof(char)*MAX_BUFFER);
    int numprocs_end = 0;

    struct timeval start_task, stop_task;
    double t_write;
    
    start_timer(start_task);
    while (1) {
      MPI_Recv(recv_buff, MAX_BUFFER, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
      if (strcmp("END", recv_buff) == 0) {
	numprocs_end++;
	if (numprocs_end == numprocs - 1) {
	  break;
	}
      } else {
	fwrite(recv_buff, strlen(recv_buff), sizeof(char), fd_out);
      }
    }
    stop_timer(start_task, stop_task, t_write);
    printf("Writer Time %0.2f\n", t_write / 1000000);
  }
    
  MPI_Barrier(MPI_COMM_WORLD);
  
  if (rank == 0) {
    stop_timer(start_w, stop_w, t_w1);
    start_timer(start_w);
  }
  
  if (rank == 0) {
    stop_timer(start, stop, t_process);
    //printf("\tEND WORKFLOW -1- (%0.2fs)\n", t_w1 / 1000000);
    //printf("\tSTART MERGE -1-\n");
    printf("======================== WORKFLOW MERGE 1 ============================\n");
    start_timer(start);
  }
  /*
  if (numprocs - 1 > 1) {
    //==========================================================================//
    //= M E R G E    M E T A E X O N - A V L     T O    2nd    W O R K F L O W =//
    //==========================================================================//
    if (rank != numprocs - 1) {
      int num_sj;        
      int numprocs_merge = numprocs - 1;
      int dec_s = 1;
      int dec_r = 2;
      int send, recv;
      int dec = 2;    
      while (1) {
	send = numprocs_merge - dec_s;
	recv = numprocs_merge - dec_r;
	if (send > 0 && recv < 0)
	  recv = 0;
	if (send <= 0 && recv <= 0)
	  break;
	for (int r = numprocs_merge - 1; r >= 0; r--) {
	  if (rank == send) {
	    //Package AVLs & Send
	    //printf("[%i]Init AVL pack send to %i\n", rank, recv);
	    int max_num_sj = 300000;
	    MPI_splice_t *mpi_avl_pack = MPI_avl_package(&max_num_sj, num_chromosomes, avls_list, &num_sj);
	    //printf("[%i]Send AVL pack to %i\n", rank, recv);
	    //sleep(1);	  
	    MPI_Send(&num_sj, 1, MPI_INT, recv, tag, MPI_COMM_WORLD);
	    MPI_Send(mpi_avl_pack, sizeof(MPI_splice_t)*num_sj, MPI_CHAR, recv, tag, MPI_COMM_WORLD);
	    //Package Metaexon & Send	  
	    MPI_metaexon_package_t *mpi_meta_pack = MPI_metaexon_package(metaexons);
	    //sleep(1);
	    MPI_Send(&mpi_meta_pack->n_metaexons, 1, MPI_INT, recv, tag, MPI_COMM_WORLD);
	    MPI_Send(mpi_meta_pack->MPI_metaexon, sizeof(MPI_metaexon_t)*mpi_meta_pack->n_metaexons, MPI_CHAR, recv, tag, MPI_COMM_WORLD);	  
	    MPI_Send(&mpi_meta_pack->n_starts_ends, 1, MPI_INT, recv, tag, MPI_COMM_WORLD);
	    MPI_Send(mpi_meta_pack->starts_ends, sizeof(MPI_breaks_t)*mpi_meta_pack->n_starts_ends, MPI_CHAR, recv, tag, MPI_COMM_WORLD);	  

	    //free(mpi_avl_pack);
	    //free(mpi_meta_pack->MPI_metaexon);
	    //free(mpi_meta_pack->starts_ends);
	  }
	  if (rank == recv) {
	    //Recv AVLs
	    int recv_num_sj;
	    //printf("[%i] recv from %i\n", rank, send);
	    MPI_Recv(&recv_num_sj, 1, MPI_INT, send, tag, MPI_COMM_WORLD, &status);
	    MPI_splice_t MPI_splice[recv_num_sj];
	    MPI_Recv(MPI_splice, sizeof(MPI_splice_t)*recv_num_sj, MPI_CHAR, send, tag, MPI_COMM_WORLD, &status);	  
	    int recv_n_metaexons;
	    MPI_Recv(&recv_n_metaexons, 1, MPI_INT, send, tag, MPI_COMM_WORLD, &status);
	    MPI_metaexon_t MPI_metaexon[recv_n_metaexons];
	    MPI_Recv(MPI_metaexon, sizeof(MPI_metaexon_t)*recv_n_metaexons, MPI_CHAR, send, tag, MPI_COMM_WORLD, &status);
	    int recv_n_starts_ends;
	    MPI_Recv(&recv_n_starts_ends, 1, MPI_INT, send, tag, MPI_COMM_WORLD, &status);
	    MPI_breaks_t MPI_breaks[recv_n_starts_ends];
	    MPI_Recv(MPI_breaks, sizeof(MPI_breaks_t)*recv_n_starts_ends, MPI_CHAR, send, tag, MPI_COMM_WORLD, &status);	  
	    //sleep(1);
	    //printf("[%i]Merge %i sp:\n", rank, recv_num_sj);
	    //Merge AVLs
	    avl_node_t *node_avl_start, *node_avl_end;
	    for (int sj = 0; sj < recv_num_sj; sj++) {
	      MPI_splice_t MPI_sj = MPI_splice[sj];
	      //printf("[%i]Merge (%i) %i: %lu - %lu\n", rank, MPI_sj.strand, MPI_sj.chr, MPI_sj.start, MPI_sj.end);
	      allocate_start_node(MPI_sj.chr - 1,
				  MPI_sj.strand,
				  MPI_sj.start,
				  MPI_sj.end,
				  MPI_sj.start,
				  MPI_sj.end,
				  FROM_READ,
				  1,
				  "",
				  &node_avl_start,
				  &node_avl_end,
				  avls_list);
	      //Refresh MPI_sj.reads_number in nodes
	      //end_data_t *end_data = node_avl_start->data;	    
	    }
	    //Merge Metaexon
	    //sleep(1);
	    //printf("[%i]Merge %i META:\n", rank, recv_n_metaexons);
	    int se_pos = 0;
	    for (int m = 0; m < recv_n_metaexons; m++) {
	      MPI_metaexon_t mpi_meta = MPI_metaexon[m];
	      //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
	      for (int s = 0; s < mpi_meta.n_starts; s++) {
		MPI_breaks_t bk = MPI_breaks[se_pos++];
		//printf("\tM-STARTS-(%i)%lu\n", bk.strand, bk.pos);
		cp_avltree *avl = avls_list->avls[bk.strand][mpi_meta.chromosome].avl;
		avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
		assert(avl_node);
		metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
				METAEXON_RIGHT_END, avl_node, metaexons);
	      }
	      for (int s = 0; s < mpi_meta.n_ends; s++) {
		MPI_breaks_t bk = MPI_breaks[se_pos++];
		//printf("\t[%i]M-ENDS-(%i)%lu\n", rank, bk.strand, bk.pos);
		cp_avltree *avl = avls_list->ends_avls[bk.strand][mpi_meta.chromosome].avl;
		avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
		assert(avl_node);
		metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
				METAEXON_LEFT_END, avl_node, metaexons);
	      }
	    }
	  }
	  send -= dec;
	  recv -= dec;
	  if (send <= 0) { send = -1; }
	}
	dec   *= 2;
	dec_s *= 2;
	dec_r *= 2;
      }
    }

    int num_sj;
    int recv_num_sj, recv_n_metaexons, recv_n_starts_ends;
    MPI_splice_t *MPI_splice;
    MPI_metaexon_package_t *mpi_meta_pack;
    MPI_metaexon_t *MPI_metaexon;
    MPI_breaks_t *MPI_breaks;

    //printf("[%i]Bcast prepare\n", rank);
    if (rank == 0) {
      int max_num_sj = 300000;
      MPI_splice  = MPI_avl_package(&max_num_sj, num_chromosomes, avls_list, &num_sj);
      mpi_meta_pack = MPI_metaexon_package(metaexons);
      recv_n_metaexons = mpi_meta_pack->n_metaexons;
      recv_n_starts_ends = mpi_meta_pack->n_starts_ends;
      MPI_metaexon = mpi_meta_pack->MPI_metaexon;
      MPI_breaks = mpi_meta_pack->starts_ends;
    }

    //printf("[%i]Bcast Send\n", rank);
    MPI_Bcast(&num_sj, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&recv_n_metaexons, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&recv_n_starts_ends, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank != 0) {
      MPI_splice  = (MPI_splice_t*)malloc(sizeof(MPI_splice_t)*num_sj);
      MPI_metaexon  = (MPI_metaexon_t *)malloc(sizeof(MPI_metaexon_t)*recv_n_metaexons);
      MPI_breaks = (MPI_breaks_t *)malloc(sizeof(MPI_breaks_t)*recv_n_starts_ends);
    }

    MPI_Bcast(MPI_splice, sizeof(MPI_splice_t)*num_sj, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(MPI_metaexon, sizeof(MPI_metaexon_t)*recv_n_metaexons, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(MPI_breaks, sizeof(MPI_breaks_t)*recv_n_starts_ends, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank != numprocs - 1) {
      //printf("[%i]Merge\n", rank);
      //Merge AVLs
      avl_node_t *node_avl_start, *node_avl_end;
      for (int sj = 0; sj < num_sj; sj++) {
	MPI_splice_t MPI_sj = MPI_splice[sj];
	//printf("[%i]Merge (%i) %i: %lu - %lu\n", rank, MPI_sj.strand, MPI_sj.chr, MPI_sj.start, MPI_sj.end);
	allocate_start_node(MPI_sj.chr - 1,
			    MPI_sj.strand,
			    MPI_sj.start,
			    MPI_sj.end,
			    MPI_sj.start,
			    MPI_sj.end,
			    FROM_READ,
			    1,
			    "",
			    &node_avl_start,
			    &node_avl_end,
			    avls_list);
	//Refresh MPI_sj.reads_number in nodes
	end_data_t *end_data = node_avl_start->data;	    
      }

      //printf("[%i]Merge AVL End\n", rank);
    
      //Merge Metaexon
      int se_pos = 0;
      for (int m = 0; m < recv_n_metaexons; m++) {
	MPI_metaexon_t mpi_meta = MPI_metaexon[m];
	//printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
	for (int s = 0; s < mpi_meta.n_starts; s++) {
	  MPI_breaks_t bk = MPI_breaks[se_pos++];
	  //printf("\tM-STARTS-(%i)%lu\n", bk.strand, bk.pos);
	  cp_avltree *avl = avls_list->avls[bk.strand][mpi_meta.chromosome].avl;
	  avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	  assert(avl_node);
	  metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			  METAEXON_RIGHT_END, avl_node, metaexons);
	}
	for (int s = 0; s < mpi_meta.n_ends; s++) {
	  MPI_breaks_t bk = MPI_breaks[se_pos++];
	  //printf("\t[%i]M-ENDS-(%i)%lu\n", rank, bk.strand, bk.pos);
	  cp_avltree *avl = avls_list->ends_avls[bk.strand][mpi_meta.chromosome].avl;
	  avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	  assert(avl_node);
	  metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			  METAEXON_LEFT_END, avl_node, metaexons);
	}
      }

      //printf("[%i]Merge Meta End\n", rank);
    
      free(MPI_splice);
      free(MPI_metaexon);
      free(MPI_breaks);
    }
  }
  */
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0) {
    stop_timer(start_w, stop_w, t_m1);
    //printf("\tSTOP MERGE -1- (%0.2fs)\n", t_m1 / 1000000);
    //printf("\tSTART WORKFLOW -2-\n");
    printf("======================== WORKFLOW 2 ============================\n");
    start_timer(start_w);
  }

  
  if (rank != numprocs - 1) {
    //====            W2 - W O R K I N G               ====//
    //==== C R E A T E    P O L L I N G    T H R E A D ====//
    MPI_data.queue_tasks = queue_sa_tasks;
    printf("batch\n");
    wf_input_file_t *wf_input_file = wf_input_file_new(NULL, batch);   
    printf("batch\n");
	
    workflow_t *wf_last = workflow_new();
    workflow_stage_function_t stage_functions_last[] = {rna_last_stage, post_pair_stage_w2_w3, MPI_output};
    char *stage_labels_last[] = {"RNA LAST STAGE", "POST PAIR", "MPI OUTPUT"};
    workflow_set_stages(3, (workflow_stage_function_t *)&stage_functions_last, stage_labels_last, wf_last);
    workflow_set_producer((workflow_producer_function_t *)file_reader, "Buffer reader", wf_last);
    workflow_set_consumer((workflow_consumer_function_t *)MPI_consumer, "SAM writer", wf_last);


    size_t num_cpus = get_optimal_cpu_num_threads();
    num_cpus -= 2;

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    cpu_set_t cpu_set;
    CPU_ZERO(&cpu_set);
    CPU_SET(num_cpus - 1, &cpu_set);
    sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set), &cpu_set);

    int ret;
    pthread_t thread;
    if ((ret = pthread_create(&thread, &attr, mpi_thread_polling, (void *)&MPI_data))) {
      printf("ERROR; return code from pthread_create() is %d\n", ret);
      exit(-1);
    }

    //=============================//
    //= W O R K E R - T H R E A D =//
    //=============================//
  
    char *file;
    int victim;
    work_package_t work;
    w1_tasks = 0;
    t_task = 0;
    
    //start_timer(start_loop);
    while (1) {
      //Process all tasks in the queue
      while (1) {
	pthread_mutex_lock(&(MPI_data.queue_mutex));
	if (MPI_data.queue_tasks->size) {
	  file = (char *)linked_list_remove_first(MPI_data.queue_tasks);
	} else {
	  file = NULL;
	}
	pthread_mutex_unlock(&(MPI_data.queue_mutex));
	if (!file) { break; }      

	start_timer(start_file);
	//====       W O R K E R    P R O C E S S    F I L E       ====//
	//=============================================================//
	char file_path[1024];
	char buffer_sa_path[1024];
	char buffer_hc_path[1024];
      
	//sprintf(file_path, "%s/buffer_sa.%s", TMP_PATH_FILES, file);
	sprintf(buffer_sa_path, "%s/buffer_sa.%s", TMP_PATH_BUFFERS, file);
	sprintf(buffer_hc_path, "%s/buffer_hc.%s", TMP_PATH_BUFFERS, file);
      
	//printf("rank.%i process file -> %s\n", rank, file_path);

	//Run workflow
	f_sa = fopen(buffer_sa_path, "r+b");
	if (f_sa == NULL) {
	  LOG_FATAL("Error opening file 'buffer_sa.tmp' \n");
	}
	sw_input.f_sa = f_sa;
      
	f_hc = fopen(buffer_hc_path, "a+b");
	if (f_hc == NULL) {
	  LOG_FATAL("Error opening file 'buffer_hc.tmp' \n");
	}
	
	wf_input_file->file = f_sa;

	//printf("[%i]Process file %s\n", rank, buffer_sa_path);
	start_timer(start_task);
	sw_input.f_hc = f_hc;
	t_file = 0;
	workflow_run_with(options->num_cpu_threads, wf_input_file, wf_last);
	//workflow_run_with(6, wf_input_file, wf_last);
	tot_t_file += t_file;
	w1_tasks++;
	stop_timer(start_task, stop_task, t_task);
	//printf("[%i]Process end\n", rank, buffer_sa_path);
	
	start_timer(start_send);
	//Send Batches to Write
	MPI_Request req;
	size_t len_send = (array_list_size(data_out->buffer) * MAX_BUFFER) / 1024; //KB
	len_send = len_send / 1024; //MB
	
	for (int i = 0; i < array_list_size(data_out->buffer); i++) {
	  batch_out_t *bo = array_list_get(i, data_out->buffer);
	  MPI_Send(bo->data, MAX_BUFFER, MPI_CHAR, numprocs - 1, 1, MPI_COMM_WORLD);
	  free_batch_out(bo);
	}
	array_list_clear(data_out->buffer, (void *)NULL);
	stop_timer(start_send, stop_send, t_send);
	
	wf_last->completed_producer = 0;
	
	fclose(f_sa);
	fclose(f_hc);
      
	//==== W O R K E R    E N D    P R O C E S S    F I L E ====//
	//==========================================================//
	stop_timer(start_file, stop_file, t_file);
	//printf("\t[%i]Task '%s' completed in %0.2f\n", rank, buffer_sa_path, t_file / 1000000);
	t_file = 0;
      }
      //Now I am idle. So, I send message to thread MPI comunicator
      work.source = rank;
    
      if (MPI_data.run) {
	work.mode = SAME_NODE;
	MPI_Send(&work, sizeof(work_package_t), MPI_CHAR, work.source, 1, MPI_COMM_WORLD);
      } else {
	work.mode = KILL_TH;
	MPI_Send(&work, sizeof(work_package_t), MPI_CHAR, work.source, 1, MPI_COMM_WORLD);
	break;
      }
      
      while (MPI_data.block) {}
      MPI_data.block = 1;
    
      //printf("\t[%i]Break block!\n", rank);
    }
    
    // free attribute and wait for the other threads
    
    void *pth_status;
    pthread_attr_destroy(&attr);
    
    if ((ret = pthread_join(thread, &pth_status))) {
      printf("ERROR; return code from pthread_join() is %d\n", ret);
      exit(-1);
    }
    
    
    //End Writer
    char *end = (char *)malloc(sizeof(char)*MAX_BUFFER);
    strcpy(end, "END\0");
    MPI_Send(end, MAX_BUFFER, MPI_CHAR, numprocs - 1, 1, MPI_COMM_WORLD);

    printf("[%i] Total Tasks %i process in %0.2fs (Time x task = %0.2fs), Total reads process %lu\n", rank, w1_tasks, t_task / 1000000, (t_task / 1000000) / w1_tasks, mpi_tot_reads);
  } else {
    //=============================//
    //= W R I T E R - T H R E A D =//
    //=============================//
    //fd_out = fopen(output_filename, "a");
    char *recv_buff = (char *)malloc(sizeof(char)*MAX_BUFFER);
    int numprocs_end = 0;

    struct timeval start_write, stop_write;
    double t_write;
    
    start_timer(start_write);
    while (1) {
      MPI_Recv(recv_buff, MAX_BUFFER, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
      if (strcmp("END", recv_buff) == 0) {
	numprocs_end++;
	if (numprocs_end == numprocs - 1) {
	  break;
	}
      } else {
	fwrite(recv_buff, strlen(recv_buff), sizeof(char), fd_out);
      }            
    }
    stop_timer(start_write, stop_write, t_write);
    printf("Writer Time %0.2f\n", t_write / 1000000);
    
  }


  MPI_Barrier(MPI_COMM_WORLD);

  
  if (rank == 0) {
    stop_timer(start_w, stop_w, t_w2);
    start_timer(start_w);
  }

  //struct timeval m_start, m_stop;
  //double m_time;
  
  if (rank == 0) {
    //printf("WORKFLOWS END Done!\n");
    stop_timer(start, stop, t_process);
    //printf("\tSTOP WORKFLOW -2- (%0.2fs), Total Task (%0.2fs)\n", t_w2 / 1000000, t_task / 1000000);
    //printf("\tSTART MERGE -2-\n");
    //printf("======================== WORKFLOW MERGE 2 ============================\n");
    start_timer(start);
  }
  /*
  if (numprocs - 1 > 1) {
    //==========================================================================//
    //= M E R G E    M E T A E X O N - A V L     T O    2nd    W O R K F L O W =//
    //==========================================================================//
    if (rank != numprocs - 1) {
      //== Package structures ==//
      int num_sj;        
      int numprocs_merge = numprocs - 1;
      int dec_s = 1;
      int dec_r = 2;
      int send, recv;
      int dec = 2;    
      while (1) {
	send = numprocs_merge - dec_s;
	recv = numprocs_merge - dec_r;
	if (send > 0 && recv < 0)
	  recv = 0;
	if (send <= 0 && recv <= 0)
	  break;
	for (int r = numprocs_merge - 1; r >= 0; r--) {
	  if (rank == send) {
	    //Package AVLs & Send
	    //printf("[%i]Init AVL pack send to %i\n", rank, recv);
	    int max_num_sj = 300000;
	    MPI_splice_t *mpi_avl_pack = MPI_avl_package(&max_num_sj, num_chromosomes, avls_list, &num_sj);
	    //printf("[%i]Send AVL pack to %i\n", rank, recv);
	    //sleep(1);	  
	    MPI_Send(&num_sj, 1, MPI_INT, recv, tag, MPI_COMM_WORLD);
	    MPI_Send(mpi_avl_pack, sizeof(MPI_splice_t)*num_sj, MPI_CHAR, recv, tag, MPI_COMM_WORLD);
	    //Package Metaexon & Send	  
	    MPI_metaexon_package_t *mpi_meta_pack = MPI_metaexon_package(metaexons);
	    //sleep(1);
	    MPI_Send(&mpi_meta_pack->n_metaexons, 1, MPI_INT, recv, tag, MPI_COMM_WORLD);
	    MPI_Send(mpi_meta_pack->MPI_metaexon, sizeof(MPI_metaexon_t)*mpi_meta_pack->n_metaexons, MPI_CHAR, recv, tag, MPI_COMM_WORLD);	  
	    MPI_Send(&mpi_meta_pack->n_starts_ends, 1, MPI_INT, recv, tag, MPI_COMM_WORLD);
	    MPI_Send(mpi_meta_pack->starts_ends, sizeof(MPI_breaks_t)*mpi_meta_pack->n_starts_ends, MPI_CHAR, recv, tag, MPI_COMM_WORLD);	  

	    //free(mpi_avl_pack);
	    //free(mpi_meta_pack->MPI_metaexon);
	    //free(mpi_meta_pack->starts_ends);
	  }
	  if (rank == recv) {
	    //Recv AVLs
	    int recv_num_sj;
	    //printf("[%i] recv from %i\n", rank, send);
	    MPI_Recv(&recv_num_sj, 1, MPI_INT, send, tag, MPI_COMM_WORLD, &status);
	    MPI_splice_t MPI_splice[recv_num_sj];
	    MPI_Recv(MPI_splice, sizeof(MPI_splice_t)*recv_num_sj, MPI_CHAR, send, tag, MPI_COMM_WORLD, &status);	  
	    int recv_n_metaexons;
	    MPI_Recv(&recv_n_metaexons, 1, MPI_INT, send, tag, MPI_COMM_WORLD, &status);
	    MPI_metaexon_t MPI_metaexon[recv_n_metaexons];
	    MPI_Recv(MPI_metaexon, sizeof(MPI_metaexon_t)*recv_n_metaexons, MPI_CHAR, send, tag, MPI_COMM_WORLD, &status);
	    int recv_n_starts_ends;
	    MPI_Recv(&recv_n_starts_ends, 1, MPI_INT, send, tag, MPI_COMM_WORLD, &status);
	    MPI_breaks_t MPI_breaks[recv_n_starts_ends];
	    MPI_Recv(MPI_breaks, sizeof(MPI_breaks_t)*recv_n_starts_ends, MPI_CHAR, send, tag, MPI_COMM_WORLD, &status);	  
	    //sleep(1);
	    //printf("[%i]Merge %i sp:\n", rank, recv_num_sj);
	    //Merge AVLs
	    avl_node_t *node_avl_start, *node_avl_end;
	    for (int sj = 0; sj < recv_num_sj; sj++) {
	      MPI_splice_t MPI_sj = MPI_splice[sj];
	      //printf("[%i]Merge (%i) %i: %lu - %lu\n", rank, MPI_sj.strand, MPI_sj.chr, MPI_sj.start, MPI_sj.end);
	      allocate_start_node(MPI_sj.chr - 1,
				  MPI_sj.strand,
				  MPI_sj.start,
				  MPI_sj.end,
				  MPI_sj.start,
				  MPI_sj.end,
				  FROM_READ,
				  1,
				  "",
				  &node_avl_start,
				  &node_avl_end,
				  avls_list);
	      //Refresh MPI_sj.reads_number in nodes
	      //end_data_t *end_data = node_avl_start->data;	    
	    }
	    //Merge Metaexon
	    //sleep(1);
	    //printf("[%i]Merge %i META:\n", rank, recv_n_metaexons);
	    int se_pos = 0;
	    for (int m = 0; m < recv_n_metaexons; m++) {
	      MPI_metaexon_t mpi_meta = MPI_metaexon[m];
	      //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
	      for (int s = 0; s < mpi_meta.n_starts; s++) {
		MPI_breaks_t bk = MPI_breaks[se_pos++];
		//printf("\tM-STARTS-(%i)%lu\n", bk.strand, bk.pos);
		cp_avltree *avl = avls_list->avls[bk.strand][mpi_meta.chromosome].avl;
		avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
		assert(avl_node);
		metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
				METAEXON_RIGHT_END, avl_node, metaexons);
	      }
	      for (int s = 0; s < mpi_meta.n_ends; s++) {
		MPI_breaks_t bk = MPI_breaks[se_pos++];
		//printf("\t[%i]M-ENDS-(%i)%lu\n", rank, bk.strand, bk.pos);
		cp_avltree *avl = avls_list->ends_avls[bk.strand][mpi_meta.chromosome].avl;
		avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
		assert(avl_node);
		metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
				METAEXON_LEFT_END, avl_node, metaexons);
	      }
	    }
	  }
	  send -= dec;
	  recv -= dec;
	  if (send <= 0) { send = -1; }
	}
	dec   *= 2;
	dec_s *= 2;
	dec_r *= 2;
      }
    }

    int num_sj;
    int recv_num_sj, recv_n_metaexons, recv_n_starts_ends;
    MPI_splice_t *MPI_splice;
    MPI_metaexon_package_t *mpi_meta_pack;
    MPI_metaexon_t *MPI_metaexon;
    MPI_breaks_t *MPI_breaks;

    
   //printf("[%i]Bcast prepare\n", rank);
    if (rank == 0) {
      int max_num_sj = 300000;
      MPI_splice  = MPI_avl_package(&max_num_sj, num_chromosomes, avls_list, &num_sj);
      mpi_meta_pack = MPI_metaexon_package(metaexons);
      recv_n_metaexons = mpi_meta_pack->n_metaexons;
      recv_n_starts_ends = mpi_meta_pack->n_starts_ends;
      MPI_metaexon = mpi_meta_pack->MPI_metaexon;
      MPI_breaks = mpi_meta_pack->starts_ends;
    }

    //printf("[%i]Bcast Send\n", rank);
    MPI_Bcast(&num_sj, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&recv_n_metaexons, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&recv_n_starts_ends, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank != 0) {
      MPI_splice  = (MPI_splice_t*)malloc(sizeof(MPI_splice_t)*num_sj);
      MPI_metaexon  = (MPI_metaexon_t *)malloc(sizeof(MPI_metaexon_t)*recv_n_metaexons);
      MPI_breaks = (MPI_breaks_t *)malloc(sizeof(MPI_breaks_t)*recv_n_starts_ends);
    }

    MPI_Bcast(MPI_splice, sizeof(MPI_splice_t)*num_sj, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(MPI_metaexon, sizeof(MPI_metaexon_t)*recv_n_metaexons, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(MPI_breaks, sizeof(MPI_breaks_t)*recv_n_starts_ends, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank != numprocs - 1) {
      //printf("[%i]Merge\n", rank);
      //Merge AVLs
      avl_node_t *node_avl_start, *node_avl_end;
      for (int sj = 0; sj < num_sj; sj++) {
	MPI_splice_t MPI_sj = MPI_splice[sj];
	//printf("[%i]Merge (%i) %i: %lu - %lu\n", rank, MPI_sj.strand, MPI_sj.chr, MPI_sj.start, MPI_sj.end);
	allocate_start_node(MPI_sj.chr - 1,
			    MPI_sj.strand,
			    MPI_sj.start,
			    MPI_sj.end,
			    MPI_sj.start,
			    MPI_sj.end,
			    FROM_READ,
			    1,
			    "",
			    &node_avl_start,
			    &node_avl_end,
			    avls_list);
	//Refresh MPI_sj.reads_number in nodes
	end_data_t *end_data = node_avl_start->data;	    
      }

      //printf("[%i]Merge AVL End\n", rank);
    
      //Merge Metaexon
      int se_pos = 0;
      for (int m = 0; m < recv_n_metaexons; m++) {
	MPI_metaexon_t mpi_meta = MPI_metaexon[m];
	//printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
	for (int s = 0; s < mpi_meta.n_starts; s++) {
	  MPI_breaks_t bk = MPI_breaks[se_pos++];
	  //printf("\tM-STARTS-(%i)%lu\n", bk.strand, bk.pos);
	  cp_avltree *avl = avls_list->avls[bk.strand][mpi_meta.chromosome].avl;
	  avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	  assert(avl_node);
	  metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			  METAEXON_RIGHT_END, avl_node, metaexons);
	}
	for (int s = 0; s < mpi_meta.n_ends; s++) {
	  MPI_breaks_t bk = MPI_breaks[se_pos++];
	  //printf("\t[%i]M-ENDS-(%i)%lu\n", rank, bk.strand, bk.pos);
	  cp_avltree *avl = avls_list->ends_avls[bk.strand][mpi_meta.chromosome].avl;
	  avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
	  assert(avl_node);
	  metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
			  METAEXON_LEFT_END, avl_node, metaexons);
	}
      }

      //printf("[%i]Merge Meta End\n", rank);
    
      free(MPI_splice);
      free(MPI_metaexon);
      free(MPI_breaks);
    }    
  }
  */
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0) {
    stop_timer(start_w, stop_w, t_m2);
    //printf("\tSTOP MERGE -2- (%0.2fs)\n", t_m2 / 1000000);
    //printf("\tSTART WORKFLOW -3-\n");
    printf("======================== WORKFLOW 3 ============================\n");
    start_timer(start_w);
  }

  
  //printf("[%i]Workflow start\n", rank);
  //== W3 ==//
  if (rank != numprocs - 1) {
    //====            W3 - W O R K I N G               ====//
    //==== C R E A T E    P O L L I N G    T H R E A D ====//
    MPI_data.queue_tasks = queue_hc_tasks;
    wf_input_file_t *wf_input_file = wf_input_file_new(NULL, batch);     

    workflow_t *wf_hc = workflow_new();
    workflow_stage_function_t stage_functions_hc[] = {rna_last_hc_stage, post_pair_stage_w2_w3, MPI_output};
    char *stage_labels_hc[] = {"RNA HARD CLIPPINGS", "POST PAIR", "MPI OUTPUT"};
    workflow_set_stages(3, (workflow_stage_function_t *)&stage_functions_hc, stage_labels_hc, wf_hc);
    workflow_set_producer((workflow_producer_function_t *)file_reader_2, "Buffer reader", wf_hc);
    workflow_set_consumer((workflow_consumer_function_t *)MPI_consumer, "SAM writer", wf_hc);

    size_t num_cpus = get_optimal_cpu_num_threads();
    num_cpus -= 2;

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    cpu_set_t cpu_set;
    CPU_ZERO(&cpu_set);
    CPU_SET(num_cpus - 1, &cpu_set);
    sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set), &cpu_set);

    int ret;
    pthread_t thread;
    if ((ret = pthread_create(&thread, &attr, mpi_thread_polling, (void *)&MPI_data))) {
      printf("ERROR; return code from pthread_create() is %d\n", ret);
      exit(-1);
    }
  
    //=============================//
    //= W O R K E R - T H R E A D =//
    //=============================//
    char *file;
    int victim;
    work_package_t work;

    //start_timer(start_loop);
    while (1) {
      //Process all tasks in the queue
      while (1) {
	pthread_mutex_lock(&(MPI_data.queue_mutex));
	if (MPI_data.queue_tasks->size) {
	  file = (char *)linked_list_remove_first(MPI_data.queue_tasks);
	} else {
	  file = NULL;
	}
	pthread_mutex_unlock(&(MPI_data.queue_mutex));
	if (!file) { break; }      

	start_timer(start_file);
	//====       W O R K E R    P R O C E S S    F I L E       ====//
	//=============================================================//
	char buffer_hc_path[1024];
      
	//sprintf(file_path, "%s/buffer_sa.%s", TMP_PATH_FILES, file);
	sprintf(buffer_hc_path, "%s/buffer_hc.%s", TMP_PATH_BUFFERS, file);
      
	//printf("rank.%i process file -> %s\n", rank, file_path);
	//Run workflow
      
	f_hc = fopen(buffer_hc_path, "r+b");
	if (f_hc == NULL) {
	  LOG_FATAL("Error opening file 'buffer_hc.tmp' \n");
	}
	
	wf_input_file->file = f_hc;

	//printf("[%i]Process file %s\n", rank, buffer_sa_path);
	sw_input.f_hc = f_hc;
	t_file = 0;
	workflow_run_with(options->num_cpu_threads, wf_input_file, wf_hc);
	tot_t_file += t_file;
	
	//Send Batches to Write
	MPI_Request req;
	size_t len_send = (array_list_size(data_out->buffer) * MAX_BUFFER) / 1024; //KB
	len_send = len_send / 1024; //MB
	
	for (int i = 0; i < array_list_size(data_out->buffer); i++) {
	  batch_out_t *bo = array_list_get(i, data_out->buffer);
	  MPI_Send(bo->data, MAX_BUFFER, MPI_CHAR, numprocs - 1, 1, MPI_COMM_WORLD);
	  free_batch_out(bo);
	}
	array_list_clear(data_out->buffer, (void *)NULL);
	
	wf_hc->completed_producer = 0;

	fclose(f_hc);
      
	//==== W O R K E R    E N D    P R O C E S S    F I L E ====//
	//==========================================================//
	stop_timer(start_file, stop_file, t_file);
	//printf("\t[%i]Task '%s' completed in %0.2f\n", rank, buffer_hc_path, t_file / 1000000);
	t_file = 0;
      }
      //Now I am idle. So, I send message to thread MPI comunicator
      work.source = rank;

      if (MPI_data.run) {
	work.mode = SAME_NODE;
	MPI_Send(&work, sizeof(work_package_t), MPI_CHAR, work.source, 1, MPI_COMM_WORLD);
      } else {
	work.mode = KILL_TH;
	MPI_Send(&work, sizeof(work_package_t), MPI_CHAR, work.source, 1, MPI_COMM_WORLD);
	break;
      }
      
      while (MPI_data.block) {}
      MPI_data.block = 1;

      //printf("\t[%i]Break block!\n", rank);
    }
    
    // free attribute and wait for the other threads
    void *pth_status;
    pthread_attr_destroy(&attr);
    
    if ((ret = pthread_join(thread, &pth_status))) {
      printf("ERROR; return code from pthread_join() is %d\n", ret);
      exit(-1);
    }
    
    //End Writer
    char *end = (char *)malloc(sizeof(char)*MAX_BUFFER);
    strcpy(end, "END\0");
    MPI_Send(end, MAX_BUFFER, MPI_CHAR, numprocs - 1, 1, MPI_COMM_WORLD);

  } else {
    //=============================//
    //= W R I T E R - T H R E A D =//
    //=============================//
    //fd_out = fopen(output_filename, "a");
    char *recv_buff = (char *)malloc(sizeof(char)*MAX_BUFFER);
    int numprocs_end = 0;
    while (1) {
      MPI_Recv(recv_buff, MAX_BUFFER, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
      if (strcmp("END", recv_buff) == 0) {
	numprocs_end++;
	if (numprocs_end == numprocs - 1) {
	  break;
	}
      } else {
	fwrite(recv_buff, strlen(recv_buff), sizeof(char), fd_out);
      }            
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0) {
    stop_timer(start_w, stop_w, t_w3);
    start_timer(start_w);
  }
  //struct timeval m_start, m_stop;
  //double m_time;
  if (rank == 0) {
    //printf("WORKFLOWS END Done!\n");
    stop_timer(start, stop, t_process);
    //printf("======================== MERGE WORKFLOW 3 ============================\n");
    start_timer(start);
  }

  //============================================================//
  //= M E R G E    M E T A E X O N - A V L    T O    W R I T E =//
  //============================================================//
  //== Package structures ==//
  if (rank != numprocs - 1) {
  /*
    if (numprocs - 1 > 1) {  
      int num_sj;        
      int numprocs_merge = numprocs - 1;
      int dec_s = 1;
      int dec_r = 2;
      int send, recv;
      int dec = 2;    
      while (1) {
	send = numprocs_merge - dec_s;
	recv = numprocs_merge - dec_r;
	if (send > 0 && recv < 0)
	  recv = 0;
	if (send <= 0 && recv <= 0)
	  break;
	for (int r = numprocs_merge - 1; r >= 0; r--) {
	  if (rank == send) {
	    //Package AVLs & Send
	    //printf("[%i]Init AVL pack send to %i\n", rank, recv);
	    int max_num_sj = 300000;
	    MPI_splice_t *mpi_avl_pack = MPI_avl_package(&max_num_sj, num_chromosomes, avls_list, &num_sj);
	    //printf("[%i]Send AVL pack to %i\n", rank, recv);
	    //sleep(1);	  
	    MPI_Send(&num_sj, 1, MPI_INT, recv, tag, MPI_COMM_WORLD);
	    MPI_Send(mpi_avl_pack, sizeof(MPI_splice_t)*num_sj, MPI_CHAR, recv, tag, MPI_COMM_WORLD);
	    //Package Metaexon & Send	  
	    MPI_metaexon_package_t *mpi_meta_pack = MPI_metaexon_package(metaexons);
	    //sleep(1);
	    MPI_Send(&mpi_meta_pack->n_metaexons, 1, MPI_INT, recv, tag, MPI_COMM_WORLD);
	    MPI_Send(mpi_meta_pack->MPI_metaexon, sizeof(MPI_metaexon_t)*mpi_meta_pack->n_metaexons, MPI_CHAR, recv, tag, MPI_COMM_WORLD);	  
	    MPI_Send(&mpi_meta_pack->n_starts_ends, 1, MPI_INT, recv, tag, MPI_COMM_WORLD);
	    MPI_Send(mpi_meta_pack->starts_ends, sizeof(MPI_breaks_t)*mpi_meta_pack->n_starts_ends, MPI_CHAR, recv, tag, MPI_COMM_WORLD);	  
	  }
	  if (rank == recv) {
	    //Recv AVLs
	    int recv_num_sj;
	    //printf("[%i] recv from %i\n", rank, send);
	    MPI_Recv(&recv_num_sj, 1, MPI_INT, send, tag, MPI_COMM_WORLD, &status);
	    MPI_splice_t MPI_splice[recv_num_sj];
	    MPI_Recv(MPI_splice, sizeof(MPI_splice_t)*recv_num_sj, MPI_CHAR, send, tag, MPI_COMM_WORLD, &status);	  
	    int recv_n_metaexons;
	    MPI_Recv(&recv_n_metaexons, 1, MPI_INT, send, tag, MPI_COMM_WORLD, &status);
	    MPI_metaexon_t MPI_metaexon[recv_n_metaexons];
	    MPI_Recv(MPI_metaexon, sizeof(MPI_metaexon_t)*recv_n_metaexons, MPI_CHAR, send, tag, MPI_COMM_WORLD, &status);
	    int recv_n_starts_ends;
	    MPI_Recv(&recv_n_starts_ends, 1, MPI_INT, send, tag, MPI_COMM_WORLD, &status);
	    MPI_breaks_t MPI_breaks[recv_n_starts_ends];
	    MPI_Recv(MPI_breaks, sizeof(MPI_breaks_t)*recv_n_starts_ends, MPI_CHAR, send, tag, MPI_COMM_WORLD, &status);	  
	    //sleep(1);
	    //printf("[%i]Merge %i sp:\n", rank, recv_num_sj);
	    //Merge AVLs
	    avl_node_t *node_avl_start, *node_avl_end;
	    for (int sj = 0; sj < recv_num_sj; sj++) {
	      MPI_splice_t MPI_sj = MPI_splice[sj];
	      //printf("[%i]Merge (%i) %i: %lu - %lu\n", rank, MPI_sj.strand, MPI_sj.chr, MPI_sj.start, MPI_sj.end);
	      allocate_start_node(MPI_sj.chr - 1,
				  MPI_sj.strand,
				  MPI_sj.start,
				  MPI_sj.end,
				  MPI_sj.start,
				  MPI_sj.end,
				  FROM_READ,
				  1,
				  "",
				  &node_avl_start,
				  &node_avl_end,
				  avls_list);
	      //Refresh MPI_sj.reads_number in nodes
	      //end_data_t *end_data = node_avl_start->data;	    
	    }
	    //Merge Metaexon
	    //sleep(1);
	    //printf("[%i]Merge %i META:\n", rank, recv_n_metaexons);
	    int se_pos = 0;
	    for (int m = 0; m < recv_n_metaexons; m++) {
	      MPI_metaexon_t mpi_meta = MPI_metaexon[m];
	      //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
	      for (int s = 0; s < mpi_meta.n_starts; s++) {
		MPI_breaks_t bk = MPI_breaks[se_pos++];
		//printf("\tM-STARTS-(%i)%lu\n", bk.strand, bk.pos);
		cp_avltree *avl = avls_list->avls[bk.strand][mpi_meta.chromosome].avl;
		avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
		assert(avl_node);
		metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
				METAEXON_RIGHT_END, avl_node, metaexons);
	      }
	      for (int s = 0; s < mpi_meta.n_ends; s++) {
		MPI_breaks_t bk = MPI_breaks[se_pos++];
		//printf("\t[%i]M-ENDS-(%i)%lu\n", rank, bk.strand, bk.pos);
		cp_avltree *avl = avls_list->ends_avls[bk.strand][mpi_meta.chromosome].avl;
		avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
		assert(avl_node);
		metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
				METAEXON_LEFT_END, avl_node, metaexons);
	      }
	    }
	  }
	  send -= dec;
	  recv -= dec;
	  if (send <= 0) { send = -1; }
	}
	dec   *= 2;
	dec_s *= 2;
	dec_r *= 2;
      }
    }
  */
    if (rank == 0) {
      double t_tmp = 0;
      struct timeval start_tmp, stop_tmp;
      start_timer(start_tmp);
      write_chromosome_avls(NULL, exact_filename, num_chromosomes, avls_list);
      stop_timer(start_tmp, stop_tmp, t_tmp);

      //printf("Chromosome Avl in %0.2f\n", t_tmp / 1000000);
    }
    
    // free memory
    //workflow_free(wf);
    //wf_input_free(wf_input);
    //Write chromosome avls
    
    // free memory
    if (bwt_index) { bwt_index_free(bwt_index); }
    metaexons_free(metaexons);
    
    batch_free(batch);
    
    //
    // end of workflow management
    //--------------------------------------------------------------------------------------
    
    if (genome)
      genome_free(genome);
    
    free(output_filename);
    free(exact_filename);
  } else {
    double t_tmp = 0;
    struct timeval start_tmp, stop_tmp;
    start_timer(start_tmp);
    fclose(fd_out);
    stop_timer(start_tmp, stop_tmp, t_tmp);
    //printf("fwrite in %0.2f\n", t_tmp / 1000000);
  }

  if (rank == 0) {
    stop_timer(start_w, stop_w, t_m3);
    //printf("\tSTOP MERGE -3 antes BARRIER- (%0.2fs)\n", t_m3 / 1000000);
    start_timer(start_w);
  }

  
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0) {
    stop_timer(start_w, stop_w, t_m3);
    //printf("\tSTOP MERGE -3- (%0.2fs)\n", t_m3 / 1000000);
  }

 exit:
  
  if (rank == 0) {
    //stop_timer(start, stop, t_merge);
    stop_timer(start_total, stop_total, t_total);
    printf("==================================================================\n");
    //printf(" MPI Genome Load Time          : %0.2f(s)\n", t_genomes / 1000000);
    //printf("      Split File Time          : %0.2f(s)\n", t_split / 1000000);

    printf(" MPI W1                        : %0.2f(s)\n", t_w1 / 1000000);
    printf(" MPI Merge-1                   : %0.2f(s)\n", t_m1 / 1000000);    
    printf(" MPI W2                        : %0.2f(s)\n", t_w2 / 1000000);
    printf(" MPI Merge-2                   : %0.2f(s)\n", t_m2 / 1000000);    
    printf(" MPI W3                        : %0.2f(s)\n", t_w3 / 1000000);
    printf(" MPI Merge-3                   : %0.2f(s)\n", t_m3 / 1000000);
    
    printf(" TOTAL TIME                    : %0.2f(s)\n", t_total / 1000000);    
    printf("==================================================================\n\n\n");
  }

  //printf("[%i]END\n", rank);
  
  MPI_Finalize();
  
}

