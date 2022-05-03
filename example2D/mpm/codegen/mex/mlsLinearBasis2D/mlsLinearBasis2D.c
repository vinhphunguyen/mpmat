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
#include "mlsLinearBasis2D_emxutil.h"
#include "inv.h"
#include "mlsLinearBasis2D_data.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 35, "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m" };

static emlrtRTEInfo emlrtRTEI = { 1, 18, "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m" };

static emlrtRTEInfo b_emlrtRTEI = { 16, 1, "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m" };

static emlrtBCInfo emlrtBCI = { -1, -1, 20, 10, "index", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

static emlrtDCInfo emlrtDCI = { 21, 15, "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 1 };

static emlrtBCInfo b_emlrtBCI = { -1, -1, 21, 15, "node", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

static emlrtBCInfo c_emlrtBCI = { -1, -1, 22, 15, "node", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

static emlrtBCInfo d_emlrtBCI = { -1, -1, 25, 31, "di", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

static emlrtDCInfo b_emlrtDCI = { 37, 14, "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 1 };

static emlrtBCInfo e_emlrtBCI = { -1, -1, 37, 14, "node", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

static emlrtBCInfo f_emlrtBCI = { -1, -1, 37, 14, "index", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

static emlrtBCInfo g_emlrtBCI = { -1, -1, 39, 3, "phi", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

static emlrtBCInfo h_emlrtBCI = { -1, -1, 39, 23, "w", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

static emlrtBCInfo i_emlrtBCI = { -1, -1, 31, 3, "w", "mlsLinearBasis2D",
  "/Users/vinhphunguyen/my-codes/MAMP/mls/mlsLinearBasis2D.m", 0 };

/* Function Definitions */
void mlsLinearBasis2D(const emlrtStack *sp, const real_T pt[2], const
                      emxArray_real_T *b_index, const emxArray_real_T *node,
                      const emxArray_real_T *di, const char_T form[2],
                      emxArray_real_T *phi)
{
  int32_T nn;
  real_T A[9];
  emxArray_real_T *w;
  int32_T c_index;
  int32_T k;
  int32_T m;
  real_T y;
  real_T wi;
  real_T b_node[2];
  real_T x[2];
  real_T scale;
  real_T absxk;
  real_T t;
  real_T r;
  real_T dv0[9];
  int32_T b_k;
  real_T p[3];
  real_T invA[9];
  int32_T b_m;
  real_T piT[3];
  real_T b_y[3];
  int32_T c_m;
  int32_T d_m;
  emlrtStack st;
  (void)form;
  st.prev = sp;
  st.tls = sp->tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

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
  nn = b_index->size[0] - 1;
  memset(&A[0], 0, 9U * sizeof(real_T));
  emxInit_real_T(sp, &w, 1, &b_emlrtRTEI, true);
  c_index = w->size[0];
  w->size[0] = b_index->size[0];
  emxEnsureCapacity(sp, (emxArray__common *)w, c_index, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  k = b_index->size[0];
  for (c_index = 0; c_index < k; c_index++) {
    w->data[c_index] = 0.0;
  }

  c_index = phi->size[0];
  phi->size[0] = b_index->size[0];
  emxEnsureCapacity(sp, (emxArray__common *)phi, c_index, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  k = b_index->size[0];
  for (c_index = 0; c_index < k; c_index++) {
    phi->data[c_index] = 0.0;
  }

  m = 0;
  while (m <= nn) {
    c_index = b_index->size[0];
    k = m + 1;
    emlrtDynamicBoundsCheckR2012b(k, 1, c_index, &emlrtBCI, sp);
    y = b_index->data[m];
    c_index = node->size[0];
    if (y == (int32_T)muDoubleScalarFloor(y)) {
      k = (int32_T)y;
    } else {
      k = (int32_T)emlrtIntegerCheckR2012b(y, &emlrtDCI, sp);
    }

    emlrtDynamicBoundsCheckR2012b(k, 1, c_index, &b_emlrtBCI, sp);
    c_index = node->size[0];
    k = (int32_T)b_index->data[m];
    emlrtDynamicBoundsCheckR2012b(k, 1, c_index, &c_emlrtBCI, sp);

    /* wi   = computeCircleSpline(pt,[xi yi],di(idm),form); */
    wi = 0.0;
    b_node[0] = node->data[(int32_T)b_index->data[m] - 1];
    b_node[1] = node->data[((int32_T)b_index->data[m] + node->size[0]) - 1];
    for (c_index = 0; c_index < 2; c_index++) {
      x[c_index] = pt[c_index] - b_node[c_index];
    }

    y = 0.0;
    scale = 2.2250738585072014E-308;
    for (k = 0; k < 2; k++) {
      absxk = muDoubleScalarAbs(x[k]);
      if (absxk > scale) {
        t = scale / absxk;
        y = 1.0 + y * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * muDoubleScalarSqrt(y);
    c_index = di->size[0];
    k = (int32_T)b_index->data[m];
    emlrtDynamicBoundsCheckR2012b(k, 1, c_index, &d_emlrtBCI, sp);
    r = y / di->data[(int32_T)b_index->data[m] - 1];
    if (r <= 1.0) {
      wi = ((1.0 - 6.0 * (r * r)) + 8.0 * muDoubleScalarPower(r, 3.0)) - 3.0 *
        muDoubleScalarPower(r, 4.0);
    }

    dv0[0] = 1.0;
    dv0[3] = node->data[(int32_T)b_index->data[m] - 1];
    dv0[6] = node->data[((int32_T)b_index->data[m] + node->size[0]) - 1];
    dv0[1] = node->data[(int32_T)b_index->data[m] - 1];
    dv0[4] = node->data[(int32_T)b_index->data[m] - 1] * node->data[(int32_T)
      b_index->data[m] - 1];
    dv0[7] = node->data[(int32_T)b_index->data[m] - 1] * node->data[((int32_T)
      b_index->data[m] + node->size[0]) - 1];
    dv0[2] = node->data[((int32_T)b_index->data[m] + node->size[0]) - 1];
    dv0[5] = node->data[(int32_T)b_index->data[m] - 1] * node->data[((int32_T)
      b_index->data[m] + node->size[0]) - 1];
    dv0[8] = node->data[((int32_T)b_index->data[m] + node->size[0]) - 1] *
      node->data[((int32_T)b_index->data[m] + node->size[0]) - 1];
    for (c_index = 0; c_index < 3; c_index++) {
      for (k = 0; k < 3; k++) {
        A[k + 3 * c_index] += wi * dv0[k + 3 * c_index];
      }
    }

    c_index = w->size[0];
    k = 1 + m;
    if ((k >= 1) && (k < c_index)) {
      b_k = k;
    } else {
      b_k = emlrtDynamicBoundsCheckR2012b(k, 1, c_index, &i_emlrtBCI, sp);
    }

    w->data[b_k - 1] = wi;
    m++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  p[0] = 1.0;
  p[1] = pt[0];
  p[2] = pt[1];
  st.site = &emlrtRSI;
  inv(&st, A, invA);
  m = 0;
  while (m <= nn) {
    c_index = b_index->size[0];
    if (m + 1 < c_index) {
      b_m = m + 1;
    } else {
      b_m = emlrtDynamicBoundsCheckR2012b(m + 1, 1, c_index, &f_emlrtBCI, sp);
    }

    y = b_index->data[b_m - 1];
    c_index = node->size[0];
    if (y == (int32_T)muDoubleScalarFloor(y)) {
      k = (int32_T)y;
    } else {
      k = (int32_T)emlrtIntegerCheckR2012b(y, &b_emlrtDCI, sp);
    }

    emlrtDynamicBoundsCheckR2012b(k, 1, c_index, &e_emlrtBCI, sp);
    k = (int32_T)b_index->data[m];
    c_index = (int32_T)b_index->data[m];
    piT[0] = 1.0;
    piT[1] = node->data[k - 1];
    piT[2] = node->data[(c_index + node->size[0]) - 1];
    y = 0.0;
    for (c_index = 0; c_index < 3; c_index++) {
      b_y[c_index] = 0.0;
      for (k = 0; k < 3; k++) {
        b_y[c_index] += p[k] * invA[k + 3 * c_index];
      }

      y += b_y[c_index] * piT[c_index];
    }

    c_index = phi->size[0];
    k = w->size[0];
    if (m + 1 < c_index) {
      c_m = m + 1;
    } else {
      c_m = emlrtDynamicBoundsCheckR2012b(m + 1, 1, c_index, &g_emlrtBCI, sp);
    }

    if (m + 1 < k) {
      d_m = m + 1;
    } else {
      d_m = emlrtDynamicBoundsCheckR2012b(m + 1, 1, k, &h_emlrtBCI, sp);
    }

    phi->data[c_m - 1] = y * w->data[d_m - 1];
    m++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  emxFree_real_T(&w);

  /* codegen -o mlsLinearBasis2D_mex mlsLinearBasis2D -args {coder.typeof(double(0), [1 2]),coder.typeof([1 2],[Inf 1]),coder.typeof(double(0), [Inf 2]),coder.typeof(double(0), [Inf 1]),coder.typeof('aa')} */
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (mlsLinearBasis2D.c) */
