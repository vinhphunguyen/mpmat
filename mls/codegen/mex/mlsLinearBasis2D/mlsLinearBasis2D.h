/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mlsLinearBasis2D.h
 *
 * Code generation for function 'mlsLinearBasis2D'
 *
 */

#ifndef __MLSLINEARBASIS2D_H__
#define __MLSLINEARBASIS2D_H__

/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blas.h"
#include "rtwtypes.h"
#include "mlsLinearBasis2D_types.h"

/* Function Declarations */
extern void mlsLinearBasis2D(const emlrtStack *sp, const real_T pt[2], const
  real_T index_data[], const int32_T index_size[1], const real_T node[2048],
  const real_T di[1024], const char_T form[12], real_T phi_data[], int32_T
  phi_size[1]);

#endif

/* End of code generation (mlsLinearBasis2D.h) */
