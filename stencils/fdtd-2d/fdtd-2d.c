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
/* fdtd-2d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "fdtd-2d.h"


/* Array initialization. */
static
void init_array (int tmax,
		 int nx,
		 int ny,
		 DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_1D(_fict_,TMAX,tmax))
{
  int i, j;

  for (i = 0; i < tmax; i++)
    _fict_[i] = (DATA_TYPE) i;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      {
	ex[i][j] = ((DATA_TYPE) i*(j+1)) / nx;
	ey[i][j] = ((DATA_TYPE) i*(j+2)) / ny;
	hz[i][j] = ((DATA_TYPE) i*(j+3)) / nx;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int nx,
		 int ny,
		 DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("ex");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, ex[i][j]);
    }
  POLYBENCH_DUMP_END("ex");
  POLYBENCH_DUMP_FINISH;

  POLYBENCH_DUMP_BEGIN("ey");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, ey[i][j]);
    }
  POLYBENCH_DUMP_END("ey");

  POLYBENCH_DUMP_BEGIN("hz");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, hz[i][j]);
    }
  POLYBENCH_DUMP_END("hz");
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_fdtd_2d(int tmax,
		    int nx,
		    int ny,
		    DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_1D(_fict_,TMAX,tmax))
{
  int t, i, j;

  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if ((_PB_NY >= 1) && (_PB_TMAX >= 1)) {
  for (t1=0;t1<=floord(2*_PB_TMAX+_PB_NY-3,32);t1++) {
    lbp=max(ceild(t1,2),ceild(32*t1-_PB_TMAX+1,32));
    ubp=min(min(floord(_PB_TMAX+_PB_NY-2,32),floord(32*t1+_PB_NY+30,64)),t1);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=t1-t2;t3<=min(min(floord(32*t1-32*t2+31*_PB_NX,32),floord(32*t1-32*t2+_PB_NX+30,32)),floord(32*t1-32*t2+30*_PB_TMAX+31*_PB_NX-30,992));t3++) {
        if ((_PB_NX >= 2) && (_PB_NY >= 2) && (t1 == 2*t2) && (t1 == 2*t3)) {
          for (t4=16*t1;t4<=min(_PB_TMAX-1,16*t1+30);t4++) {
            if (t1%2 == 0) {
              ey[0][0] = _fict_[t4];;
            }
            for (t6=t4+1;t6<=min(16*t1+31,t4+_PB_NX-1);t6++) {
              if (t1%2 == 0) {
                ey[(-t4+t6)][0] = ey[(-t4+t6)][0] - SCALAR_VAL(0.5)*(hz[(-t4+t6)][0]-hz[(-t4+t6)-1][0]);;
              }
            }
            for (t5=t4+1;t5<=min(16*t1+31,t4+_PB_NY-1);t5++) {
              if (t1%2 == 0) {
                ex[0][(-t4+t5)] = ex[0][(-t4+t5)] - SCALAR_VAL(0.5)*(hz[0][(-t4+t5)]-hz[0][(-t4+t5)-1]);;
              }
              if (t1%2 == 0) {
                ey[0][(-t4+t5)] = _fict_[t4];;
              }
              for (t6=t4+1;t6<=min(16*t1+31,t4+_PB_NX-1);t6++) {
                if (t1%2 == 0) {
                  ey[(-t4+t6)][(-t4+t5)] = ey[(-t4+t6)][(-t4+t5)] - SCALAR_VAL(0.5)*(hz[(-t4+t6)][(-t4+t5)]-hz[(-t4+t6)-1][(-t4+t5)]);;
                }
                if (t1%2 == 0) {
                  ex[(-t4+t6)][(-t4+t5)] = ex[(-t4+t6)][(-t4+t5)] - SCALAR_VAL(0.5)*(hz[(-t4+t6)][(-t4+t5)]-hz[(-t4+t6)][(-t4+t5)-1]);;
                }
                if (t1%2 == 0) {
                  hz[(-t4+t6-1)][(-t4+t5-1)] = hz[(-t4+t6-1)][(-t4+t5-1)] - SCALAR_VAL(0.7)* (ex[(-t4+t6-1)][(-t4+t5-1)+1] - ex[(-t4+t6-1)][(-t4+t5-1)] + ey[(-t4+t6-1)+1][(-t4+t5-1)] - ey[(-t4+t6-1)][(-t4+t5-1)]);;
                }
              }
            }
          }
        }
        if ((_PB_NX >= 2) && (t1 == t2+t3) && (t1 <= 2*t2-1)) {
          for (t4=max(32*t1-32*t2,32*t2-_PB_NY+1);t4<=min(_PB_TMAX-1,32*t1-32*t2+30);t4++) {
            for (t5=32*t2;t5<=min(32*t2+31,t4+_PB_NY-1);t5++) {
              ex[0][(-t4+t5)] = ex[0][(-t4+t5)] - SCALAR_VAL(0.5)*(hz[0][(-t4+t5)]-hz[0][(-t4+t5)-1]);;
              ey[0][(-t4+t5)] = _fict_[t4];;
              for (t6=t4+1;t6<=min(32*t1-32*t2+31,t4+_PB_NX-1);t6++) {
                ey[(-t4+t6)][(-t4+t5)] = ey[(-t4+t6)][(-t4+t5)] - SCALAR_VAL(0.5)*(hz[(-t4+t6)][(-t4+t5)]-hz[(-t4+t6)-1][(-t4+t5)]);;
                ex[(-t4+t6)][(-t4+t5)] = ex[(-t4+t6)][(-t4+t5)] - SCALAR_VAL(0.5)*(hz[(-t4+t6)][(-t4+t5)]-hz[(-t4+t6)][(-t4+t5)-1]);;
                hz[(-t4+t6-1)][(-t4+t5-1)] = hz[(-t4+t6-1)][(-t4+t5-1)] - SCALAR_VAL(0.7)* (ex[(-t4+t6-1)][(-t4+t5-1)+1] - ex[(-t4+t6-1)][(-t4+t5-1)] + ey[(-t4+t6-1)+1][(-t4+t5-1)] - ey[(-t4+t6-1)][(-t4+t5-1)]);;
              }
            }
          }
        }
        if ((_PB_NX >= 2) && (_PB_NY == 1) && (t1 == 2*t2) && (t1 == 2*t3)) {
          for (t4=16*t1;t4<=min(_PB_TMAX-1,16*t1+30);t4++) {
            if (t1%2 == 0) {
              ey[0][0] = _fict_[t4];;
            }
            for (t6=t4+1;t6<=min(16*t1+31,t4+_PB_NX-1);t6++) {
              if (t1%2 == 0) {
                ey[(-t4+t6)][0] = ey[(-t4+t6)][0] - SCALAR_VAL(0.5)*(hz[(-t4+t6)][0]-hz[(-t4+t6)-1][0]);;
              }
            }
          }
        }
        if ((_PB_NX == 1) && (_PB_NY >= 2) && (t1 == 2*t2) && (t1 == 2*t3)) {
          for (t4=16*t1;t4<=min(_PB_TMAX-1,16*t1+30);t4++) {
            if (t1%2 == 0) {
              ey[0][0] = _fict_[t4];;
            }
            for (t5=t4+1;t5<=min(16*t1+31,t4+_PB_NY-1);t5++) {
              if (t1%2 == 0) {
                ex[0][(-t4+t5)] = ex[0][(-t4+t5)] - SCALAR_VAL(0.5)*(hz[0][(-t4+t5)]-hz[0][(-t4+t5)-1]);;
              }
              if (t1%2 == 0) {
                ey[0][(-t4+t5)] = _fict_[t4];;
              }
            }
          }
        }
        if ((_PB_NX == 1) && (t1 == t2+t3) && (t1 <= 2*t2-1)) {
          for (t4=max(32*t1-32*t2,32*t2-_PB_NY+1);t4<=min(_PB_TMAX-1,32*t1-32*t2+31);t4++) {
            for (t5=32*t2;t5<=min(32*t2+31,t4+_PB_NY-1);t5++) {
              ex[0][(-t4+t5)] = ex[0][(-t4+t5)] - SCALAR_VAL(0.5)*(hz[0][(-t4+t5)]-hz[0][(-t4+t5)-1]);;
              ey[0][(-t4+t5)] = _fict_[t4];;
            }
          }
        }
        if ((_PB_NX <= 1) && (_PB_NY == 1) && (t1 == 2*t2) && (t1 == 2*t3)) {
          for (t4=16*t1;t4<=min(_PB_TMAX-1,16*t1+31);t4++) {
            if (t1%2 == 0) {
              ey[0][0] = _fict_[t4];;
            }
          }
        }
        if ((_PB_NX == 0) && (_PB_NY >= 2) && (t1 == t2+t3)) {
          for (t4=max(32*t1-32*t2,32*t2-_PB_NY+1);t4<=min(_PB_TMAX-1,32*t1-32*t2+31);t4++) {
            for (t5=max(32*t2,t4);t5<=min(32*t2+31,t4+_PB_NY-1);t5++) {
              ey[0][(-t4+t5)] = _fict_[t4];;
            }
          }
        }
        if ((_PB_NX == 1) && (_PB_NY >= 2) && (t1 == 2*t2) && (t1 == 2*t3) && (t1 <= floord(_PB_TMAX-32,16))) {
          if (t1%2 == 0) {
            ey[0][0] = _fict_[(16*t1+31)];;
          }
        }
        if ((_PB_NX >= 2) && (t1 == t2+t3) && (t1 <= min(floord(32*t2+_PB_TMAX-32,32),2*t2-1))) {
          for (t5=32*t2;t5<=min(32*t2+31,32*t1-32*t2+_PB_NY+30);t5++) {
            ex[0][(-32*t1+32*t2+t5-31)] = ex[0][(-32*t1+32*t2+t5-31)] - SCALAR_VAL(0.5)*(hz[0][(-32*t1+32*t2+t5-31)]-hz[0][(-32*t1+32*t2+t5-31)-1]);;
            ey[0][(-32*t1+32*t2+t5-31)] = _fict_[(32*t1-32*t2+31)];;
          }
        }
        if ((_PB_NX >= 2) && (t1 == 2*t2) && (t1 == 2*t3) && (t1 <= floord(_PB_TMAX-32,16))) {
          if (t1%2 == 0) {
            ey[0][0] = _fict_[(16*t1+31)];;
          }
        }
        if ((_PB_NY >= 2) && (t1 == 2*t2) && (t1 <= 2*t3-2)) {
          for (t4=max(16*t1,32*t3-_PB_NX+1);t4<=min(_PB_TMAX-1,16*t1+30);t4++) {
            for (t6=32*t3;t6<=min(32*t3+31,t4+_PB_NX-1);t6++) {
              if (t1%2 == 0) {
                ey[(-t4+t6)][0] = ey[(-t4+t6)][0] - SCALAR_VAL(0.5)*(hz[(-t4+t6)][0]-hz[(-t4+t6)-1][0]);;
              }
            }
            for (t5=t4+1;t5<=min(16*t1+31,t4+_PB_NY-1);t5++) {
              for (t6=32*t3;t6<=min(32*t3+31,t4+_PB_NX-1);t6++) {
                if (t1%2 == 0) {
                  ey[(-t4+t6)][(-t4+t5)] = ey[(-t4+t6)][(-t4+t5)] - SCALAR_VAL(0.5)*(hz[(-t4+t6)][(-t4+t5)]-hz[(-t4+t6)-1][(-t4+t5)]);;
                }
                if (t1%2 == 0) {
                  ex[(-t4+t6)][(-t4+t5)] = ex[(-t4+t6)][(-t4+t5)] - SCALAR_VAL(0.5)*(hz[(-t4+t6)][(-t4+t5)]-hz[(-t4+t6)][(-t4+t5)-1]);;
                }
                if (t1%2 == 0) {
                  hz[(-t4+t6-1)][(-t4+t5-1)] = hz[(-t4+t6-1)][(-t4+t5-1)] - SCALAR_VAL(0.7)* (ex[(-t4+t6-1)][(-t4+t5-1)+1] - ex[(-t4+t6-1)][(-t4+t5-1)] + ey[(-t4+t6-1)+1][(-t4+t5-1)] - ey[(-t4+t6-1)][(-t4+t5-1)]);;
                }
              }
            }
          }
        }
        if (t1 <= min(2*t2-1,t2+t3-1)) {
          for (t4=max(max(32*t1-32*t2,32*t2-_PB_NY+1),32*t3-_PB_NX+1);t4<=min(_PB_TMAX-1,32*t1-32*t2+31);t4++) {
            for (t5=32*t2;t5<=min(32*t2+31,t4+_PB_NY-1);t5++) {
              for (t6=32*t3;t6<=min(32*t3+31,t4+_PB_NX-1);t6++) {
                ey[(-t4+t6)][(-t4+t5)] = ey[(-t4+t6)][(-t4+t5)] - SCALAR_VAL(0.5)*(hz[(-t4+t6)][(-t4+t5)]-hz[(-t4+t6)-1][(-t4+t5)]);;
                ex[(-t4+t6)][(-t4+t5)] = ex[(-t4+t6)][(-t4+t5)] - SCALAR_VAL(0.5)*(hz[(-t4+t6)][(-t4+t5)]-hz[(-t4+t6)][(-t4+t5)-1]);;
                hz[(-t4+t6-1)][(-t4+t5-1)] = hz[(-t4+t6-1)][(-t4+t5-1)] - SCALAR_VAL(0.7)* (ex[(-t4+t6-1)][(-t4+t5-1)+1] - ex[(-t4+t6-1)][(-t4+t5-1)] + ey[(-t4+t6-1)+1][(-t4+t5-1)] - ey[(-t4+t6-1)][(-t4+t5-1)]);;
              }
            }
          }
        }
        if ((_PB_NY == 1) && (t1 == 2*t2) && (t1 <= 2*t3-2)) {
          for (t4=max(16*t1,32*t3-_PB_NX+1);t4<=min(_PB_TMAX-1,16*t1+31);t4++) {
            for (t6=32*t3;t6<=min(32*t3+31,t4+_PB_NX-1);t6++) {
              if (t1%2 == 0) {
                ey[(-t4+t6)][0] = ey[(-t4+t6)][0] - SCALAR_VAL(0.5)*(hz[(-t4+t6)][0]-hz[(-t4+t6)-1][0]);;
              }
            }
          }
        }
        if ((_PB_NY >= 2) && (t1 == 2*t2) && (t1 <= min(floord(_PB_TMAX-32,16),2*t3-2))) {
          for (t6=32*t3;t6<=min(32*t3+31,16*t1+_PB_NX+30);t6++) {
            if (t1%2 == 0) {
              ey[(-16*t1+t6-31)][0] = ey[(-16*t1+t6-31)][0] - SCALAR_VAL(0.5)*(hz[(-16*t1+t6-31)][0]-hz[(-16*t1+t6-31)-1][0]);;
            }
          }
        }
      }
      if (_PB_NX <= -1) {
        for (t4=max(32*t1-32*t2,32*t2-_PB_NY+1);t4<=min(_PB_TMAX-1,32*t1-32*t2+31);t4++) {
          for (t5=max(32*t2,t4);t5<=min(32*t2+31,t4+_PB_NY-1);t5++) {
            ey[0][(-t4+t5)] = _fict_[t4];;
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
  int tmax = TMAX;
  int nx = NX;
  int ny = NY;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(ex,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_2D_ARRAY_DECL(ey,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_2D_ARRAY_DECL(hz,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_1D_ARRAY_DECL(_fict_,DATA_TYPE,TMAX,tmax);

  /* Initialize array(s). */
  init_array (tmax, nx, ny,
	      POLYBENCH_ARRAY(ex),
	      POLYBENCH_ARRAY(ey),
	      POLYBENCH_ARRAY(hz),
	      POLYBENCH_ARRAY(_fict_));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  _mm_mfence();
  kernel_fdtd_2d (tmax, nx, ny,
		  POLYBENCH_ARRAY(ex),
		  POLYBENCH_ARRAY(ey),
		  POLYBENCH_ARRAY(hz),
		  POLYBENCH_ARRAY(_fict_));


  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(nx, ny, POLYBENCH_ARRAY(ex),
				    POLYBENCH_ARRAY(ey),
				    POLYBENCH_ARRAY(hz)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(ex);
  POLYBENCH_FREE_ARRAY(ey);
  POLYBENCH_FREE_ARRAY(hz);
  POLYBENCH_FREE_ARRAY(_fict_);

  return 0;
}
