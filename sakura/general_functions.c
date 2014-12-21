#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "utils.h"


void* sakura_malloc(size_t items, size_t size, char *message){

  void *ptr = malloc(items*size);

  if(ptr == 0){
    printf("Out of memory at %s\n", message);
  }

  return ptr;
}


void* sakura_calloc(size_t items, size_t size, char* message){

  void *ptr = calloc(items, size);

  if(ptr == 0){
    printf("Out of memory %s\n", message);
  }

  return ptr;
}
