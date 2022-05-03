/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * inv.h
 *
 * Code generation for function 'inv'
 *
 */

#ifndef __INV_H__
#define __INV_H__

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
extern void inv(const emlrtStack *sp, const real_T x[9], real_T y[9]);

#endif

/* End of code generation (inv.h) */
