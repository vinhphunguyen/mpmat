/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mlsLinearBasis2D.c
 *
 * Code generation for function 'mlsLinearBasis2D'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mlsLinearBasis2D.h"
#include "eml_warning.h"
#include "mlsLinearBasis2D_data.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 23, "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m" };

static emlrtRSInfo b_emlrtRSI = { 31, "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m" };

static emlrtRSInfo c_emlrtRSI = { 16, "computeCircleSpline",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/computeCircleSpline.m" };

static emlrtRSInfo d_emlrtRSI = { 1, "mldivide",
  "/Applications/MATLAB_R2015a.app/toolbox/eml/lib/matlab/ops/mldivide.p" };

static emlrtRSInfo e_emlrtRSI = { 48, "eml_lusolve",
  "/Applications/MATLAB_R2015a.app/toolbox/eml/lib/matlab/eml/eml_lusolve.m" };

static emlrtRSInfo f_emlrtRSI = { 235, "eml_lusolve",
  "/Applications/MATLAB_R2015a.app/toolbox/eml/lib/matlab/eml/eml_lusolve.m" };

static emlrtRSInfo g_emlrtRSI = { 76, "eml_lusolve",
  "/Applications/MATLAB_R2015a.app/toolbox/eml/lib/matlab/eml/eml_lusolve.m" };

static emlrtMCInfo emlrtMCI = { 20, 5, "error",
  "/Applications/MATLAB_R2015a.app/toolbox/eml/lib/matlab/lang/error.m" };

static emlrtBCInfo emlrtBCI = { -1, -1, 20, 10, "index", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

static emlrtDCInfo emlrtDCI = { 21, 15, "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 1 };

static emlrtBCInfo b_emlrtBCI = { 1, 1024, 21, 15, "node", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

static emlrtDCInfo b_emlrtDCI = { 35, 14, "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 1 };

static emlrtBCInfo c_emlrtBCI = { 1, 1024, 35, 14, "node", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

static emlrtBCInfo d_emlrtBCI = { -1, -1, 35, 14, "index", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

static emlrtBCInfo e_emlrtBCI = { -1, -1, 37, 3, "phi", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

static emlrtBCInfo f_emlrtBCI = { -1, -1, 37, 19, "w", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

static emlrtBCInfo g_emlrtBCI = { -1, -1, 26, 3, "w", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

static emlrtRSInfo h_emlrtRSI = { 20, "error",
  "/Applications/MATLAB_R2015a.app/toolbox/eml/lib/matlab/lang/error.m" };

/* Function Declarations */
static void error(const emlrtStack *sp, const mxArray *b, emlrtMCInfo *location);

/* Function Definitions */
static void error(const emlrtStack *sp, const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(sp, 0, NULL, 1, &pArray, "error", true, location);
}

void mlsLinearBasis2D(const emlrtStack *sp, const real_T pt[2], const real_T
                      index_data[], const int32_T index_size[1], const real_T
                      node[2048], const real_T di[1024], const char_T form[12],
                      real_T phi_data[], int32_T phi_size[1])
{
  int32_T nn;
  real_T A[9];
  int32_T w_size_idx_0;
  int32_T r1;
  int32_T r2;
  real_T w_data[52];
  int32_T r3;
  real_T b_node[2];
  real_T x[2];
  real_T maxval;
  real_T a21;
  real_T absxk;
  real_T t;
  real_T r;
  boolean_T b_bool;
  int32_T exitg1;
  static const char_T cv0[12] = { 'c', 'u', 'b', 'i', 'c', '_', 's', 'p', 'l',
    'i', 'n', 'e' };

  int32_T c_bool;
  real_T wi;
  real_T dv0[9];
  int32_T b_r2;
  real_T p[3];
  int32_T rtemp;
  real_T c[3];
  int32_T b_r3;
  int32_T c_r3;
  int32_T d_r3;
  const mxArray *y;
  char_T u[28];
  const mxArray *m0;
  static const int32_T iv0[2] = { 1, 28 };

  static const char_T varargin_1[28] = { 'G', 'r', 'r', '.', ' ', 'U', 'n', 'k',
    'n', 'o', 'w', 'n', ' ', 'f', 'u', 'n', 'c', 't', 'i', 'o', 'n', 'a', 'l',
    ' ', 'f', 'o', 'r', 'm' };

  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;

  /*  Compute the MLS shape function at point pt for all nodes within the */
  /*  support of this point pt. */
  /*  Basis used is linear basis pT = [1 x y] */
  /*  */
  /*  Vinh Phu Nguyen */
  /*  The University of Adelaide, Australia */
  /*  9 October 2015. */
  /*  -------------------------------------- */
  /*       compute the moment matrix A */
  /*  -------------------------------------- */
  nn = index_size[0] - 1;
  memset(&A[0], 0, 9U * sizeof(real_T));
  w_size_idx_0 = index_size[0];
  r1 = index_size[0];
  for (r2 = 0; r2 < r1; r2++) {
    w_data[r2] = 0.0;
  }

  phi_size[0] = index_size[0];
  r1 = index_size[0];
  for (r2 = 0; r2 < r1; r2++) {
    phi_data[r2] = 0.0;
  }

  r3 = 0;
  while (r3 <= nn) {
    r2 = r3 + 1;
    emlrtDynamicBoundsCheckR2012b(r2, 1, index_size[0], &emlrtBCI, sp);
    if (index_data[r3] == (int32_T)muDoubleScalarFloor(index_data[r3])) {
      r2 = (int32_T)index_data[r3];
    } else {
      r2 = (int32_T)emlrtIntegerCheckR2012b(index_data[r3], &emlrtDCI, sp);
    }

    emlrtDynamicBoundsCheckR2012b(r2, 1, 1024, &b_emlrtBCI, sp);
    st.site = &emlrtRSI;

    /*  Compute cubic and quartic spline function */
    /*  Inputs: */
    /*  x (1x2)  : coordinate of point at which w is to be evaluated */
    /*  xI (1x2) : coord of node I */
    /*  d        : size of the support */
    b_node[0] = node[(int32_T)index_data[r3] - 1];
    b_node[1] = node[(int32_T)index_data[r3] + 1023];
    for (r2 = 0; r2 < 2; r2++) {
      x[r2] = pt[r2] - b_node[r2];
    }

    maxval = 0.0;
    a21 = 2.2250738585072014E-308;
    for (r1 = 0; r1 < 2; r1++) {
      absxk = muDoubleScalarAbs(x[r1]);
      if (absxk > a21) {
        t = a21 / absxk;
        maxval = 1.0 + maxval * t * t;
        a21 = absxk;
      } else {
        t = absxk / a21;
        maxval += t * t;
      }
    }

    maxval = a21 * muDoubleScalarSqrt(maxval);
    r = maxval / di[(int32_T)index_data[r3] - 1];
    b_bool = false;
    r1 = 0;
    do {
      exitg1 = 0;
      if (r1 < 12) {
        if (form[r1] != cv0[r1]) {
          exitg1 = 1;
        } else {
          r1++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);

    if (b_bool) {
      c_bool = 0;
    } else {
      c_bool = -1;
    }

    switch (c_bool) {
     case 0:
      /*  Compute cubic spline function */
      if (r <= 0.5) {
        maxval = (0.66666666666666663 - 4.0 * r * r) + 4.0 * r * r * r;
      } else if ((r > 0.5) && (r <= 1.0)) {
        maxval = ((1.3333333333333333 - 4.0 * r) + 4.0 * r * r) -
          1.3333333333333333 * r * r * r;
      } else {
        maxval = 0.0;
      }

      wi = maxval;
      break;

     case 1:
      /*  Compute cubic spline function */
      if (r <= 1.0) {
        wi = ((1.0 - 6.0 * r * r) + 8.0 * r * r * r) - 3.0 * r * r * r * r;
      } else {
        wi = 0.0;
      }
      break;

     default:
      b_st.site = &c_emlrtRSI;
      for (r2 = 0; r2 < 28; r2++) {
        u[r2] = varargin_1[r2];
      }

      y = NULL;
      m0 = emlrtCreateCharArray(2, iv0);
      emlrtInitCharArrayR2013a(&b_st, 28, m0, &u[0]);
      emlrtAssign(&y, m0);
      c_st.site = &h_emlrtRSI;
      error(&c_st, y, &emlrtMCI);
      break;
    }

    dv0[0] = 1.0;
    dv0[3] = node[(int32_T)index_data[r3] - 1];
    dv0[6] = node[(int32_T)index_data[r3] + 1023];
    dv0[1] = node[(int32_T)index_data[r3] - 1];
    dv0[4] = node[(int32_T)index_data[r3] - 1] * node[(int32_T)index_data[r3] -
      1];
    dv0[7] = node[(int32_T)index_data[r3] - 1] * node[(int32_T)index_data[r3] +
      1023];
    dv0[2] = node[(int32_T)index_data[r3] + 1023];
    dv0[5] = node[(int32_T)index_data[r3] - 1] * node[(int32_T)index_data[r3] +
      1023];
    dv0[8] = node[(int32_T)index_data[r3] + 1023] * node[(int32_T)index_data[r3]
      + 1023];
    for (r2 = 0; r2 < 3; r2++) {
      for (r1 = 0; r1 < 3; r1++) {
        A[r1 + 3 * r2] += wi * dv0[r1 + 3 * r2];
      }
    }

    r2 = 1 + r3;
    if ((r2 >= 1) && (r2 < w_size_idx_0)) {
      b_r2 = r2;
    } else {
      b_r2 = emlrtDynamicBoundsCheckR2012b(r2, 1, w_size_idx_0, &g_emlrtBCI, sp);
    }

    w_data[b_r2 - 1] = wi;
    r3++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  p[0] = 1.0;
  p[1] = pt[0];
  p[2] = pt[1];
  st.site = &b_emlrtRSI;
  b_st.site = &d_emlrtRSI;
  c_st.site = &e_emlrtRSI;
  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = muDoubleScalarAbs(A[0]);
  a21 = muDoubleScalarAbs(A[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if (muDoubleScalarAbs(A[2]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }

  A[r2] /= A[r1];
  A[r3] /= A[r1];
  A[3 + r2] -= A[r2] * A[3 + r1];
  A[3 + r3] -= A[r3] * A[3 + r1];
  A[6 + r2] -= A[r2] * A[6 + r1];
  A[6 + r3] -= A[r3] * A[6 + r1];
  if (muDoubleScalarAbs(A[3 + r3]) > muDoubleScalarAbs(A[3 + r2])) {
    rtemp = r2;
    r2 = r3;
    r3 = rtemp;
  }

  A[3 + r3] /= A[3 + r2];
  A[6 + r3] -= A[3 + r3] * A[6 + r2];
  if ((A[r1] == 0.0) || (A[3 + r2] == 0.0) || (A[6 + r3] == 0.0)) {
    d_st.site = &f_emlrtRSI;
    e_st.site = &g_emlrtRSI;
    eml_warning(&e_st);
  }

  c[1] = p[r2] - p[r1] * A[r2];
  c[2] = (p[r3] - p[r1] * A[r3]) - c[1] * A[3 + r3];
  c[2] /= A[6 + r3];
  c[0] = p[r1] - c[2] * A[6 + r1];
  c[1] -= c[2] * A[6 + r2];
  c[1] /= A[3 + r2];
  c[0] -= c[1] * A[3 + r1];
  c[0] /= A[r1];
  r3 = 0;
  while (r3 <= nn) {
    if (r3 + 1 < index_size[0]) {
      b_r3 = r3 + 1;
    } else {
      b_r3 = emlrtDynamicBoundsCheckR2012b(r3 + 1, 1, index_size[0], &d_emlrtBCI,
        sp);
    }

    maxval = index_data[b_r3 - 1];
    if (maxval == (int32_T)muDoubleScalarFloor(maxval)) {
      r2 = (int32_T)maxval;
    } else {
      r2 = (int32_T)emlrtIntegerCheckR2012b(maxval, &b_emlrtDCI, sp);
    }

    emlrtDynamicBoundsCheckR2012b(r2, 1, 1024, &c_emlrtBCI, sp);
    p[0] = 1.0;
    p[1] = node[(int32_T)index_data[r3] - 1];
    p[2] = node[(int32_T)index_data[r3] + 1023];
    maxval = 0.0;
    for (r1 = 0; r1 < 3; r1++) {
      maxval += c[r1] * p[r1];
    }

    if (r3 + 1 < phi_size[0]) {
      c_r3 = r3 + 1;
    } else {
      c_r3 = emlrtDynamicBoundsCheckR2012b(r3 + 1, 1, phi_size[0], &e_emlrtBCI,
        sp);
    }

    if (r3 + 1 < w_size_idx_0) {
      d_r3 = r3 + 1;
    } else {
      d_r3 = emlrtDynamicBoundsCheckR2012b(r3 + 1, 1, w_size_idx_0, &f_emlrtBCI,
        sp);
    }

    phi_data[c_r3 - 1] = maxval * w_data[d_r3 - 1];
    r3++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }
}

/* End of code generation (mlsLinearBasis2D.c) */
