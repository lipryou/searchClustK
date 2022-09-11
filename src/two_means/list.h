#ifndef INCLUDE_list_h_
#define INCLUDE_list_h_

#include "base.h"

typedef struct item {
  int ind;
  struct item *next;
} item;
typedef item* list;

typedef enum {false, true} boolean;

extern void pushTo(list que, int k);
extern int getTopFrom(list que);
extern boolean getkFrom(list que, int k);
extern int getDim(list que);
extern boolean isIn(list que, int k);
extern void free_que(list que);
extern list make_dummy();
extern void print_que(list que);
extern int list2vec(int *vec, list l);

#endif //INCLUDE_list_h_
