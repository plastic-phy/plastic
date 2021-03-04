#include "silly-test.h"
#include <stdio.h>
#include <stdlib.h>

int add(int fst, int snd) {
  return fst + snd;
}

sabbia* make_sabbia(void) {
  printf("lol");
  sabbia* out = (sabbia*) malloc(sizeof(sabbia));
  int pollo[] = {1, 2, 3}; 
  out->el = 69;
  return out;
}
