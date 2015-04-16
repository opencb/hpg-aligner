#include <multi/multi_aligner.h>

#define TMP_PATH     "TMP_FILES"
#define OUTPUT_PATH  "MPI_output"

int split_input_file(char *fq_str, char *tmp_path, int numprocs) {
  FILE *tmp_fd;
  FILE *fq_file = fopen(fq_str, "r");
  if (!fq_file) { printf("ERROR opening file\n"); exit(-1); }

  //Prepare PATH tmp files
  char cmd[2048];
  sprintf(cmd, "rm -rf %s", tmp_path);
  system(cmd);  
  create_directory(tmp_path);

  //Calculate the number of reads per file
  int reads_split;
  size_t reads_positions[numprocs + 1];

  fseek(fq_file, 0L, SEEK_END);
  size_t fd_total_bytes = ftell(fq_file);
  fseek(fq_file, 0L, SEEK_SET);
    
  int i;
  size_t seek_bytes = 0;
  size_t inc_bytes = fd_total_bytes / numprocs;
  size_t offset;
  char line[2048];
    
  for (i = 0; i < numprocs; i++) {
    offset = 0;
    fseek(fq_file, seek_bytes, SEEK_SET);
    fgets(line, 2048, fq_file);
    while (line[0] != '@') {
      offset += strlen(line);
      fgets(line, 2048, fq_file);
    }      
    seek_bytes += offset;
    reads_positions[i] = seek_bytes;
    seek_bytes += inc_bytes;
  }
  reads_positions[i] = fd_total_bytes;
  fclose(fq_file);

  fq_file = fopen(fq_str, "r");
  char sequence[MAX_READ_ID_LENGTH];
  size_t total_size = 0;
  char tmp_file[1024];
  char *res;

  //printf("SPLIT FILES\n");
  for (int rank = 0; rank < numprocs; rank++) {
    sprintf(tmp_file, "%s/%i.tmp", tmp_path, rank);
    tmp_fd = fopen(tmp_file, "w");
    //printf("[%i] %s => %lu < %lu\n", rank, tmp_file, total_size, reads_positions[rank + 1]);

    while (total_size < reads_positions[rank + 1]) {
      fgets(sequence, MAX_READ_ID_LENGTH, fq_file);
      fputs(sequence, tmp_fd);
      total_size += strlen(sequence);
    }
    fclose(tmp_fd);
  }
  fclose(fq_file);

}

int hpg_multialigner_main(options_t *options, int argc, char *argv[]) {

  int numprocs;
  int rank;  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &name_len);

  printf("[%i/%i] %s\n", rank, numprocs, processor_name);

  MPI_Barrier(MPI_COMM_WORLD);

  char *mapper_cli = strdup(options->command);
  char *file_input = strdup(options->in_filename);

  char pwd[2048];
  if (getcwd(pwd, sizeof(pwd)) == NULL) {
    perror("getcwd() error");
  }

  char *path_output = (char *)malloc(sizeof(char)*(strlen(OUTPUT_PATH) + strlen(pwd) + 1024));
  if (!options->output_name) {
    sprintf(path_output, "%s/%s/", pwd, OUTPUT_PATH);
  } else {
    path_output = strdup(options->output_name);  
  }

  if (rank == 0) {
    printf("=====================================================\n"); 
    printf("START Multi Aligner\n"); 
    printf("=====================================================\n"); 
    printf("MAPPER COMMAND LINE:  %s\n", mapper_cli);
    printf("INPUT  FILE        :  %s\n", file_input);
    printf("OUTPUT PATH        :  %s\n", path_output);
    printf("=====================================================\n"); 
  }
  
  //================ First split input File ====================//
  int nfiles;  
  int ntasks, extra_tasks;

  if (rank == 0) {
    printf("[%i] Spliting input file...\n", rank);
    nfiles = split_input_file(file_input, TMP_PATH, numprocs);      
    printf("[%i] Spliting END!\n", rank);
  }
  
  MPI_Bcast(&nfiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  ntasks = nfiles / numprocs;
  extra_tasks = nfiles % numprocs;
  
  if (rank < extra_tasks) {
    ntasks++;
  }
  //============================================================//


  //================= SPLIT INPUT CLI AND TRANSFORM IT ==========================//
  size_t max_len_cli = strlen(mapper_cli) + strlen(path_output) + strlen(TMP_PATH) + 1024;
  char *cli_tmp = (char *)malloc(sizeof(char)*max_len_cli);

  char *cli_in  = (char *)malloc(sizeof(char)*max_len_cli);
  char *cli_out = (char *)malloc(sizeof(char)*max_len_cli);
  char *cli_end = (char *)malloc(sizeof(char)*max_len_cli);
  
  size_t len_cli = 0;
  unsigned char find_mark = 0, first_in = 0, first_out = 0;
  
  //if (rank == 0) {
  char *str = mapper_cli;
  char *pch;
  
  pch = strtok (str," ");
  while (pch != NULL) {
    if (strcmp(pch, "%I") == 0 || strcmp(pch, "%i") == 0) {	
      strcpy(cli_in, cli_tmp);
      len_cli = 0;
      cli_tmp[0] = '\0';
      find_mark++;
      first_in = !first_out;
    } else if (strcmp(pch, "%O") == 0 || strcmp(pch, "%o") == 0) {
      strcpy(cli_out, cli_tmp);
      len_cli = 0;
      cli_tmp[0] = '\0';
      find_mark++;
      first_out = !first_in;
    } else {
      strcpy(cli_tmp + len_cli, pch);
      len_cli = strlen(cli_tmp);
      strcpy(cli_tmp + len_cli, " \0");
      len_cli++;
    }
    
    pch = strtok(NULL, " ");
    
  }
  
  strcpy(cli_end, cli_tmp);
  
  if (find_mark != 2) {
    if (rank == 0) {
      printf("ERROR in CLI Mapper\n");
      exit(-1);
    }  
  }

  //=============================================================================//


  //================================ MPI PROCESS =======================================//
  //COMMON OUTPUT PATH 
  if (rank == 0) {
    char cmd[strlen(path_output) + 1024];
    sprintf(cmd, "rm -rf %s", path_output);
    system(cmd);
    create_directory(path_output);  
    for (int i = 0; i < numprocs; i++) {
      char path_out_tmp[strlen(path_output) + 1024];
      sprintf(path_out_tmp, "%s/%i.out", path_output, i);
      printf("Create path %s\n", path_out_tmp);
      create_directory(path_out_tmp);      
    }
  }
    
  //PROCESS FILES
  MPI_Barrier(MPI_COMM_WORLD);
  
  char *mapper_run = (char *)malloc(sizeof(char)*max_len_cli);
  char node_out[strlen(path_output) + 1024];
  sprintf(node_out, "%s/%i.out", path_output, rank);

  if (first_in) {
    sprintf(mapper_run, "%s %s/%i.tmp %s %s %s", cli_in, TMP_PATH, rank, cli_out, node_out, cli_end);
  } else {
    sprintf(mapper_run, "%s %s %s %s/%i.tmp %s", cli_out, node_out, cli_in, TMP_PATH, rank, cli_end);
  }
  
  printf("[%i] : %s\n", rank, mapper_run);
  //system(...);  
  //====================================================================================//



  //================= MPI for input aligner ====================//
  //printf("[%i] Run Tasks %i\n", rank, ntasks);
  //============================================================//  


  MPI_Finalize();
  
  //printf("[%i] Finished\n", rank);

}
