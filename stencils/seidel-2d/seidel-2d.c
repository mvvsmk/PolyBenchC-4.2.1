#include <omp.h>
#include <math.h>
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* seidel-2d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "seidel-2d.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A[i][j] = ((DATA_TYPE) i*(j+2) + 2) / n;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
    }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_seidel_2d(int tsteps,
		      int n,
		      DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int t, i, j;

  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if ((_PB_N >= 3) && (_PB_TSTEPS >= 1)) {
  for (t1=0;t1<=floord(2*_PB_TSTEPS+_PB_N-4,32);t1++) {
    lbp=max(ceild(t1,2),ceild(32*t1-_PB_TSTEPS+1,32));
    ubp=min(min(floord(_PB_TSTEPS+_PB_N-3,32),floord(32*t1+_PB_N+29,64)),t1);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=max(ceild(64*t2-_PB_N-28,32),t1);t3<=min(min(min(min(floord(_PB_TSTEPS+_PB_N-3,16),floord(32*t1-32*t2+_PB_N+29,16)),floord(32*t1+_PB_N+60,32)),floord(64*t2+_PB_N+59,32)),floord(32*t2+_PB_TSTEPS+_PB_N+28,32));t3++) {
        for (t4=max(max(max(32*t1-32*t2,32*t2-_PB_N+2),16*t3-_PB_N+2),-32*t2+32*t3-_PB_N-29);t4<=min(min(min(min(_PB_TSTEPS-1,32*t2+30),16*t3+14),32*t1-32*t2+31),-32*t2+32*t3+30);t4++) {
          for (t5=max(max(32*t2,t4+1),32*t3-t4-_PB_N+2);t5<=min(min(32*t2+31,32*t3-t4+30),t4+_PB_N-2);t5++) {
            for (t6=max(32*t3,t4+t5+1);t6<=min(32*t3+31,t4+t5+_PB_N-2);t6++) {
              A[(-t4+t5)][(-t4-t5+t6)] = (A[(-t4+t5)-1][(-t4-t5+t6)-1] + A[(-t4+t5)-1][(-t4-t5+t6)] + A[(-t4+t5)-1][(-t4-t5+t6)+1] + A[(-t4+t5)][(-t4-t5+t6)-1] + A[(-t4+t5)][(-t4-t5+t6)] + A[(-t4+t5)][(-t4-t5+t6)+1] + A[(-t4+t5)+1][(-t4-t5+t6)-1] + A[(-t4+t5)+1][(-t4-t5+t6)] + A[(-t4+t5)+1][(-t4-t5+t6)+1])/SCALAR_VAL(9.0);;
            }
          }
        }
      }
    }
  }
}

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int tsteps = TSTEPS;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  _mm_mfence();
  kernel_seidel_2d (tsteps, n, POLYBENCH_ARRAY(A));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);

  return 0;
}
