/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eml_warning.c
 *
 * Code generation for function 'eml_warning'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mlsLinearBasis2D.h"
#include "eml_warning.h"

/* Variable Definitions */
static emlrtMCInfo b_emlrtMCI = { 16, 13, "eml_warning",
  "/Applications/MATLAB_R2015a.app/toolbox/eml/lib/matlab/eml/eml_warning.m" };

static emlrtMCInfo c_emlrtMCI = { 16, 5, "eml_warning",
  "/Applications/MATLAB_R2015a.app/toolbox/eml/lib/matlab/eml/eml_warning.m" };

static emlrtRSInfo i_emlrtRSI = { 16, "eml_warning",
  "/Applications/MATLAB_R2015a.app/toolbox/eml/lib/matlab/eml/eml_warning.m" };

/* Function Declarations */
static const mxArray *message(const emlrtStack *sp, const mxArray *b,
  emlrtMCInfo *location);
static void warning(const emlrtStack *sp, const mxArray *b, emlrtMCInfo
                    *location);

/* Function Definitions */
static const mxArray *message(const emlrtStack *sp, const mxArray *b,
  emlrtMCInfo *location)
{
  const mxArray *pArray;
  const mxArray *m3;
  pArray = b;
  return emlrtCallMATLABR2012b(sp, 1, &m3, 1, &pArray, "message", true, location);
}

static void warning(const emlrtStack *sp, const mxArray *b, emlrtMCInfo
                    *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(sp, 0, NULL, 1, &pArray, "warning", true, location);
}

void eml_warning(const emlrtStack *sp)
{
  int32_T i0;
  static const char_T varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T u[27];
  const mxArray *y;
  static const int32_T iv1[2] = { 1, 27 };

  const mxArray *m1;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  for (i0 = 0; i0 < 27; i0++) {
    u[i0] = varargin_1[i0];
  }

  y = NULL;
  m1 = emlrtCreateCharArray(2, iv1);
  emlrtInitCharArrayR2013a(sp, 27, m1, &u[0]);
  emlrtAssign(&y, m1);
  st.site = &i_emlrtRSI;
  warning(&st, message(&st, y, &b_emlrtMCI), &c_emlrtMCI);
}

/* End of code generation (eml_warning.c) */
