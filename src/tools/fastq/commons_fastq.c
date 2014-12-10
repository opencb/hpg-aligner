/*
 * commons_fastq.c
 *
 *  Created on: May 22, 2013
 *      Author: jtarraga
 *              
 */

#include "commons_fastq.h"

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
//--------------------------------------------------------------------
