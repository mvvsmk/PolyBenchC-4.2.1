/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
#ifndef _DOITGEN_H
# define _DOITGEN_H

/* Default to LARGE_DATASET. */
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET) && !defined(XL1_1) && !defined(XL1_2) && !defined(XL1_3) && !defined(XL1_5) && !defined(DOUBLE_XL)
#  define LARGE_DATASET
# endif

# if !defined(NQ) && !defined(NR) && !defined(NP)
/* Define sample dataset sizes. */
#  ifdef MINI_DATASET
#   define NQ 8
#   define NR 10
#   define NP 12
#  endif

#  ifdef SMALL_DATASET
#   define NQ 20
#   define NR 25
#   define NP 30
#  endif

#  ifdef MEDIUM_DATASET
#   define NQ 40
#   define NR 50
#   define NP 60
#  endif

#  ifdef LARGE_DATASET
#   define NQ 140
#   define NR 150
#   define NP 160
#  endif

#  ifdef EXTRALARGE_DATASET
#   define NQ 220
#   define NR 250
#   define NP 270
#  endif

#  ifdef XL1_1
#   define NQ 242
#   define NR 275
#   define NP 297
#  endif

#  ifdef XL1_2
#   define NQ 264
#   define NR 300
#   define NP 324
#  endif

#  ifdef XL1_3
#   define NQ 286
#   define NR 325
#   define NP 351
#  endif

#  ifdef XL1_5
#   define NQ 330
#   define NR 375
#   define NP 405
#  endif

#  ifdef DOUBLE_XL
#   define NQ 440
#   define NR 500
#   define NP 540
#  endif

#endif /* !(NQ NR NP) */

# define _PB_NQ POLYBENCH_LOOP_BOUND(NQ,nq)
# define _PB_NR POLYBENCH_LOOP_BOUND(NR,nr)
# define _PB_NP POLYBENCH_LOOP_BOUND(NP,np)


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

#endif /* !_DOITGEN_H */
