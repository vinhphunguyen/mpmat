/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * inv.c
 *
 * Code generation for function 'inv'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mlsLinearBasis2D.h"
#include "inv.h"
#include "eml_warning.h"
#include "norm.h"

/* Variable Definitions */
static emlrtRSInfo b_emlrtRSI = { 27, "inv",
  "/Applications/MATLAB_R2015a.app/toolbox/eml/lib/matlab/matfun/inv.m" };

static emlrtRSInfo c_emlrtRSI = { 40, "inv",
  "/Applications/MATLAB_R2015a.app/toolbox/eml/lib/matlab/matfun/inv.m" };

static emlrtRSInfo d_emlrtRSI = { 44, "inv",
  "/Applications/MATLAB_R2015a.app/toolbox/eml/lib/matlab/matfun/inv.m" };

static emlrtMCInfo c_emlrtMCI = { 29, 23, "eml_flt2str",
  "/Applications/MATLAB_R2015a.app/toolbox/eml/lib/matlab/eml/eml_flt2str.m" };

static emlrtMCInfo d_emlrtMCI = { 29, 15, "eml_flt2str",
  "/Applications/MATLAB_R2015a.app/toolbox/eml/lib/matlab/eml/eml_flt2str.m" };

static emlrtRSInfo f_emlrtRSI = { 29, "eml_flt2str",
  "/Applications/MATLAB_R2015a.app/toolbox/eml/lib/matlab/eml/eml_flt2str.m" };

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[14]);
static const mxArray *b_sprintf(const emlrtStack *sp, const mxArray *b, const
  mxArray *c, const mxArray *d, emlrtMCInfo *location);
static const mxArray *c_sprintf(const emlrtStack *sp, const mxArray *b, const
  mxArray *c, emlrtMCInfo *location);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *d_sprintf,
  const char_T *identifier, char_T y[14]);
static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[14]);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[14])
{
  k_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *b_sprintf(const emlrtStack *sp, const mxArray *b, const
  mxArray *c, const mxArray *d, emlrtMCInfo *location)
{
  const mxArray *pArrays[3];
  const mxArray *m5;
  pArrays[0] = b;
  pArrays[1] = c;
  pArrays[2] = d;
  return emlrtCallMATLABR2012b(sp, 1, &m5, 3, pArrays, "sprintf", true, location);
}

static const mxArray *c_sprintf(const emlrtStack *sp, const mxArray *b, const
  mxArray *c, emlrtMCInfo *location)
{
  const mxArray *pArrays[2];
  const mxArray *m6;
  pArrays[0] = b;
  pArrays[1] = c;
  return emlrtCallMATLABR2012b(sp, 1, &m6, 2, pArrays, "sprintf", true, location);
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *d_sprintf,
  const char_T *identifier, char_T y[14])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  b_emlrt_marshallIn(sp, emlrtAlias(d_sprintf), &thisId, y);
  emlrtDestroyArray(&d_sprintf);
}

static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[14])
{
  int32_T iv5[2];
  int32_T i2;
  for (i2 = 0; i2 < 2; i2++) {
    iv5[i2] = 1 + 13 * i2;
  }

  emlrtCheckBuiltInR2012b(sp, msgId, src, "char", false, 2U, iv5);
  emlrtImportCharArrayR2014b(sp, src, ret, 14);
  emlrtDestroyArray(&src);
}

void inv(const emlrtStack *sp, const real_T x[9], real_T y[9])
{
  real_T b_x[9];
  int32_T p1;
  int32_T p2;
  int32_T p3;
  real_T absx11;
  real_T absx21;
  real_T absx31;
  int32_T itmp;
  real_T b_y;
  static const char_T cv0[8] = { '%', '%', '%', 'd', '.', '%', 'd', 'e' };

  char_T u[8];
  const mxArray *c_y;
  static const int32_T iv0[2] = { 1, 8 };

  const mxArray *m0;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *f_y;
  char_T cv1[14];
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  memcpy(&b_x[0], &x[0], 9U * sizeof(real_T));
  p1 = 0;
  p2 = 3;
  p3 = 6;
  absx11 = muDoubleScalarAbs(x[0]);
  absx21 = muDoubleScalarAbs(x[1]);
  absx31 = muDoubleScalarAbs(x[2]);
  if ((absx21 > absx11) && (absx21 > absx31)) {
    p1 = 3;
    p2 = 0;
    b_x[0] = x[1];
    b_x[1] = x[0];
    b_x[3] = x[4];
    b_x[4] = x[3];
    b_x[6] = x[7];
    b_x[7] = x[6];
  } else {
    if (absx31 > absx11) {
      p1 = 6;
      p3 = 0;
      b_x[0] = x[2];
      b_x[2] = x[0];
      b_x[3] = x[5];
      b_x[5] = x[3];
      b_x[6] = x[8];
      b_x[8] = x[6];
    }
  }

  absx21 = b_x[1] / b_x[0];
  b_x[1] /= b_x[0];
  absx11 = b_x[2] / b_x[0];
  b_x[2] /= b_x[0];
  b_x[4] -= absx21 * b_x[3];
  b_x[5] -= absx11 * b_x[3];
  b_x[7] -= absx21 * b_x[6];
  b_x[8] -= absx11 * b_x[6];
  if (muDoubleScalarAbs(b_x[5]) > muDoubleScalarAbs(b_x[4])) {
    itmp = p2;
    p2 = p3;
    p3 = itmp;
    b_x[1] = absx11;
    b_x[2] = absx21;
    absx11 = b_x[4];
    b_x[4] = b_x[5];
    b_x[5] = absx11;
    absx11 = b_x[7];
    b_x[7] = b_x[8];
    b_x[8] = absx11;
  }

  absx31 = b_x[5];
  b_y = b_x[4];
  absx21 = b_x[5] / b_x[4];
  b_x[8] -= absx21 * b_x[7];
  absx11 = (absx21 * b_x[1] - b_x[2]) / b_x[8];
  absx21 = -(b_x[1] + b_x[7] * absx11) / b_x[4];
  y[p1] = ((1.0 - b_x[3] * absx21) - b_x[6] * absx11) / b_x[0];
  y[p1 + 1] = absx21;
  y[p1 + 2] = absx11;
  absx11 = -(absx31 / b_y) / b_x[8];
  absx21 = (1.0 - b_x[7] * absx11) / b_x[4];
  y[p2] = -(b_x[3] * absx21 + b_x[6] * absx11) / b_x[0];
  y[p2 + 1] = absx21;
  y[p2 + 2] = absx11;
  absx11 = 1.0 / b_x[8];
  absx21 = -b_x[7] * absx11 / b_x[4];
  y[p3] = -(b_x[3] * absx21 + b_x[6] * absx11) / b_x[0];
  y[p3 + 1] = absx21;
  y[p3 + 2] = absx11;
  st.site = &b_emlrtRSI;
  absx11 = norm(x);
  absx21 = norm(y);
  absx31 = 1.0 / (absx11 * absx21);
  if ((absx11 == 0.0) || (absx21 == 0.0) || (absx31 == 0.0)) {
    b_st.site = &c_emlrtRSI;
    eml_warning(&b_st);
  } else {
    if (muDoubleScalarIsNaN(absx31) || (absx31 < 2.2204460492503131E-16)) {
      b_st.site = &d_emlrtRSI;
      for (p1 = 0; p1 < 8; p1++) {
        u[p1] = cv0[p1];
      }

      c_y = NULL;
      m0 = emlrtCreateCharArray(2, iv0);
      emlrtInitCharArrayR2013a(&b_st, 8, m0, &u[0]);
      emlrtAssign(&c_y, m0);
      d_y = NULL;
      m0 = emlrtCreateDoubleScalar(14.0);
      emlrtAssign(&d_y, m0);
      e_y = NULL;
      m0 = emlrtCreateDoubleScalar(6.0);
      emlrtAssign(&e_y, m0);
      f_y = NULL;
      m0 = emlrtCreateDoubleScalar(absx31);
      emlrtAssign(&f_y, m0);
      c_st.site = &f_emlrtRSI;
      emlrt_marshallIn(&c_st, c_sprintf(&c_st, b_sprintf(&c_st, c_y, d_y, e_y,
        &c_emlrtMCI), f_y, &d_emlrtMCI), "sprintf", cv1);
      b_st.site = &d_emlrtRSI;
      b_eml_warning(&b_st, cv1);
    }
  }
}

/* End of code generation (inv.c) */
