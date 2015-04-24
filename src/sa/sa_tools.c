#include "sa_tools.h"

//--------------------------------------------------------------------------------------

#define PROGRESS 100000000

size_t PREFIX_TABLE_K_VALUE = 16;
size_t PREFIX_TABLE_NT_VALUE[256];

char *global_sa_genome;

//--------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------

int ncmp_by_id(void const *a, void const *b) { 
  return (strncmp(&global_sa_genome[(*(uint *)a)], &global_sa_genome[(*(uint *)b)], 1000));
}

void compute_sa(uint *sa, uint num_sufixes) {
  qsort(sa, num_sufixes, sizeof(uint), ncmp_by_id);
}

//--------------------------------------------------------------------------------------

void compute_lcp(char *s, uint *sa, uint *lcp, uint num_sufixes) {
  uint count = 0;

  char *curr, *prev;
  int h = 0;
  lcp[0] = 0;
  uint progress = 0;
  for (int i = 1; i < num_sufixes; i++) {
    progress++;
    if (progress % PROGRESS == 0) printf("\tlcp processing %0.2f %c (max. 255 = %u)...\n", 
					 100.0f * progress / num_sufixes, '%', count); 
    
    prev = s + sa[i - 1];
    curr = s + sa[i];
    h = 0;
    while (*prev == *curr) {
      h++;
      prev++;
      curr++;
    }
    lcp[i] = h;
    if (h >= 255) count++;
  }
  printf("\t, num lcp >= 255 -> %u\n", count);
}

//--------------------------------------------------------------------------------------

void compute_child(uint *lcp, int *child, uint num_sufixes) {
  const int no_value = -1; //len + 10;
  for (uint i = 0; i < num_sufixes; i++){
    child[i] = no_value;
  }

  int last_index = no_value;
  int stack_size = num_sufixes / 2;

  int *stack = (int *) malloc(stack_size * sizeof(int));
  int index = 0;

  printf("\tcomputing up and down values...(%u)\n", num_sufixes);
  // compute up and down values
  stack[index] = 0;
  for (uint i = 1; i < num_sufixes; i++) {
    //    printf("i = %i -> lcp[i] = %i, lcp[stack[index]] = %i, index = %i\n", i, lcp[i], lcp[stack[index]], index);
    //    if (i == 322) exit(-1);

    while (lcp[i] < lcp[stack[index]]) {
      last_index = stack[index];
      if (--index < 0) index = 0;
      if (lcp[i] <= lcp[stack[index]] && lcp[stack[index]] != lcp[last_index]) {
	child[stack[index]] = last_index;
	//if (i - 1 == 20603) { printf("000\n"); exit(-1); }
      }
    }
    //now LCP[i] >= LCP[top] holds
    if (last_index != no_value){
      //if (i - 1 == 20603) { printf("111\n"); exit(-1); }
      child[i-1] = last_index;
      last_index = no_value;
    }
    if (index + 1 >= stack_size) { printf("Overflow (%i)\n", index); exit(-1); }
    stack[++index] = i;
  }

  //last row (fix for last character of sequence not being unique
  //  while (0 == lcp[stack[index]] || no_value == lcp[stack[index]] ){
  while (0 < lcp[stack[index]] ){
    last_index = stack[index];
    if (--index < 0) index = 0;
    if (0 <= lcp[stack[index]] && lcp[stack[index]] != lcp[last_index]) {
      child[stack[index]] = last_index;
    }
  }
  printf("\tcomputing up and down values...Done\n");

  printf("\tcompuing next l-index values...\n");
  // compute next l-index values
  index = 0;
  stack[index] = 0;
  for (uint i = 1; i < num_sufixes; i++) {
    //printf("i = %i -> lcp[i] = %i, lcp[stack[index]] = %i, index = %i\n", i, lcp[i], lcp[stack[index]], index);
    while (lcp[i] < lcp[stack[index]]) {
      if (--index < 0) index = 0;
    }
    last_index = stack[index];
    if (lcp[i] == lcp[last_index]) {
      if (--index < 0) index = 0;
      child[last_index] = i;
    }
    if (index + 1 >= stack_size) { printf("Overflow (%i)\n", index); exit(-1); }
    stack[++index] = i;
  }
  printf("\tcompuing next l-index values...Done\n");

  free(stack);
}

//--------------------------------------------------------------------------------------

inline size_t compute_prefix_value(char *prefix, int len) {
  size_t value = 0;
  for (int i = len - 1, j = 0; i >= 0; i--, j += 2) {
    value += ((1LL * PREFIX_TABLE_NT_VALUE[(unsigned char)prefix[i]]) << j);
  }
  return value;
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

