#ifndef SILLY_TEST_H
#define SILLY_TEST_H

typedef struct s {
  int el;
  int* arr;
} sabbia;

sabbia* make_sabbia(void);
int add(int fst, int snd);

#endif
