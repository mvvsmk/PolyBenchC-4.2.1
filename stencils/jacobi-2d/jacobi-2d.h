/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
#ifndef _JACOBI_2D_H
# define _JACOBI_2D_H

/* Default to LARGE_DATASET. */
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET) && !defined(XL1_1) && !defined(XL1_2) && !defined(XL1_3) && !defined(XL1_5) && !defined(DOUBLE_XL)
#  define LARGE_DATASET
# endif

# if !defined(TSTEPS) && !defined(N)
/* Define sample dataset sizes. */
#  ifdef MINI_DATASET
#   define TSTEPS 20
#   define N 30
#  endif

#  ifdef SMALL_DATASET
#   define TSTEPS 40
#   define N 90
#  endif

#  ifdef MEDIUM_DATASET
#   define TSTEPS 100
#   define N 250
#  endif

#  ifdef LARGE_DATASET
#   define TSTEPS 500
#   define N 1300
#  endif

#  ifdef EXTRALARGE_DATASET
#   define TSTEPS 1000
#   define N 2800
#  endif

#  ifdef XL1_1
#   define TSTEPS 1100
#   define N 3080
#  endif

#  ifdef XL1_2
#   define TSTEPS 1200
#   define N 3360
#  endif

#  ifdef XL1_3
#   define TSTEPS 1300
#   define N 3640
#  endif

#  ifdef XL1_5
#   define TSTEPS 1500
#   define N 4200
#  endif

#  ifdef DOUBLE_XL
#   define TSTEPS 2000
#   define N 5600
#  endif


#endif /* !(TSTEPS N) */

# define _PB_TSTEPS POLYBENCH_LOOP_BOUND(TSTEPS,tsteps)
# define _PB_N POLYBENCH_LOOP_BOUND(N,n)


/* Default data type */
# if !defined(DATA_TYPE_IS_INT) && !defined(DATA_TYPE_IS_FLOAT) && !defined(DATA_TYPE_IS_DOUBLE)
#  define DATA_TYPE_IS_DOUBLE
# endif

#ifdef DATA_TYPE_IS_INT
#  define DATA_TYPE int
#  define DATA_PRINTF_MODIFIER "%d "
#endif

#ifdef DATA_TYPE_IS_FLOAT
#  define DATA_TYPE float
#  define DATA_PRINTF_MODIFIER "%0.2f "
#  define SCALAR_VAL(x) x##f
#  define SQRT_FUN(x) sqrtf(x)
#  define EXP_FUN(x) expf(x)
#  define POW_FUN(x,y) powf(x,y)
# endif

#ifdef DATA_TYPE_IS_DOUBLE
#  define DATA_TYPE double
#  define DATA_PRINTF_MODIFIER "%0.2lf "
#  define SCALAR_VAL(x) x
#  define SQRT_FUN(x) sqrt(x)
#  define EXP_FUN(x) exp(x)
#  define POW_FUN(x,y) pow(x,y)
# endif

#endif /* !_JACOBI_2D_H */
