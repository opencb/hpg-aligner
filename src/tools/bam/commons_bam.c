/*
 * commons_bam.c
 *
 *  Created on: May 22, 2013
 *      Author: jtarraga
 *              
 */

#include "commons_bam.h"
//------------------------------------------------------------------------

int read_progress = 0;

//--------------------------------------------------------------------

void free_argtable(int num_options, void **argtable) {
  if (argtable != NULL) {
    arg_freetable(argtable, num_options);
    free(argtable);
  }
}

//--------------------------------------------------------------------

void usage_argtable(char *exec_name, char *command_name, void **argtable) {
  printf("Usage: %s %s\n", exec_name, command_name);
  arg_print_syntaxv(stdout, argtable, "\n");
  arg_print_glossary(stdout, argtable, "%-50s\t%s\n");
  exit(-1);
}

//--------------------------------------------------------------------

int parse_range(int *min, int *max, char *range, char *msg) {

  *min = NO_VALUE;
  *max = NO_VALUE;

  if (!range) {
    return 1;
  }

  int lmin, lmax;
  char *p, *lrange = strdup(range);

  p = strstr(lrange, ",");
  if (p) {
    //lmax = atoi(p + 1);
    if (strlen(p+1) == 0) {
      lmax = NO_VALUE;
    } else {
      if (sscanf((p+1), "%d", &lmax) != 1) {
	printf("\nError: Invalid maximum value in the %s (%s)\n", msg, range);
	free(lrange);
	return 0;
      }
    }
    if (lrange == p) {
      lmin = NO_VALUE;
    } else {
      *p = '\0';
      //lmin = atoi(lrange);
      if (sscanf(lrange, "%d", &lmin) != 1) {
	printf("\nError: Invalid minimum value in the %s (%s)\n", msg, range);
	free(lrange);
	return 0;
      }
    }
  } else {
    //lmin = atoi(lrange);
    if (sscanf(lrange, "%d", &lmin) != 1) {
      printf("\nError: Invalid minimum value in the %s (%s)\n", msg, range);
      free(lrange);
      return 0;
    }
    lmax = NO_VALUE;
  }
  
  if (lmin != NO_VALUE && lmin < 0) {
    printf("\nError: Invalid %s (%s). Minimum value (%i) must be greater than 0\n",
	   msg, range, lmin);
    free(lrange);
    return 0;
  }

  if (lmax != NO_VALUE && lmax < 0) {
    printf("\nError: Invalid %s (%s). Maximum value (%i) must be greater than 0\n",
	   msg, range, lmax);
    free(lrange);
    return 0;
  }

  if (lmin != NO_VALUE && lmax != NO_VALUE && lmin > lmax) {
    printf("\nError: Invalid %s (%s). Maximum value (%i) must be greater than minimum value (%i)\n",
	   msg, range, lmax, lmin);
    free(lrange);
    return 0;
  }

  *min = lmin;
  *max = lmax;

  // free local range
  free(lrange);
  return 1;
}

//--------------------------------------------------------------------

region_table_t *build_region_table(char *bam_filename, char *by_string, 
				   char *by_gff_file) {  
  char **chromosomes;
  int num_chromosomes = 0;

  // read chromosome info from the bam file
  bam_file_t *bam_file = bam_fopen(bam_filename);
  num_chromosomes = bam_file->bam_header_p->n_targets;
  if (num_chromosomes <= 0) {
    bam_fclose(bam_file);
    return NULL;
  }

  chromosomes = (char **) calloc(num_chromosomes, sizeof(char *));
  for (int i = 0; i < num_chromosomes; i++) {
    chromosomes[i] = strdup(bam_file->bam_header_p->target_name[i]);
    //printf("%i of %i: %s\n", i, num_chromosomes, chromosomes[i]);
  }
  
  bam_fclose(bam_file);

  // create the region table from chromosome info
  region_table_t *region_table = new_region_table(num_chromosomes, chromosomes);

  if (region_table) {
    int ret = 0;
    // insert regions from string or gff file
    if (by_string) {
      ret = region_table_parse_from_string(by_string, region_table);
    } else if (by_gff_file) {
      ret = region_table_parse_from_gff_file(by_gff_file, region_table);
    }
    
    if (!ret) {
      free_region_table(region_table);
      region_table = NULL;
    } else {
      finish_region_table_loading(region_table);
    }
  }

  // return the table
  return region_table;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

