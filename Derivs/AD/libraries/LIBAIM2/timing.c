#include <time.h>
#include <unistd.h>

void gettime_(double *forttime) {
  //  *forttime = (( (double)clock() )*100 / ((double)CLOCKS_PER_SEC));
  *forttime =  (double)(clock()/CLOCKS_PER_SEC);
  sleep(1);
}

