/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_mlsLinearBasis2D_api.c
 *
 * Code generation for function '_coder_mlsLinearBasis2D_api'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mlsLinearBasis2D.h"
#include "_coder_mlsLinearBasis2D_api.h"
#include "mlsLinearBasis2D_emxutil.h"
#include "mlsLinearBasis2D_data.h"

/* Variable Definitions */
static emlrtRTEInfo c_emlrtRTEI = { 1, 1, "_coder_mlsLinearBasis2D_api", "" };

/* Function Declarations */
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *pt,
  const char_T *identifier))[2];
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[2];
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_index,
  const char_T *identifier, emxArray_real_T *y);
static const mxArray *emlrt_marshallOut(const emxArray_real_T *u);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *node, const
  char_T *identifier, emxArray_real_T *y);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *form, const
  char_T *identifier, char_T y[2]);
static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[2]);
static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[2];
static void m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[2]);

/* Function Definitions */
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *pt,
  const char_T *identifier))[2]
{
  real_T (*y)[2];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = d_emlrt_marshallIn(sp, emlrtAlias(pt), &thisId);
  emlrtDestroyArray(&pt);
  return y;
}
  static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[2]
{
  real_T (*y)[2];
  y = l_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_index,
  const char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  f_emlrt_marshallIn(sp, emlrtAlias(b_index), &thisId, y);
  emlrtDestroyArray(&b_index);
}

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *y;
  static const int32_T iv4[1] = { 0 };

  const mxArray *m3;
  y = NULL;
  m3 = emlrtCreateNumericArray(1, iv4, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m3, (void *)u->data);
  emlrtSetDimensions((mxArray *)m3, u->size, 1);
  emlrtAssign(&y, m3);
  return y;
}

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  m_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *node, const
  char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  h_emlrt_marshallIn(sp, emlrtAlias(node), &thisId, y);
  emlrtDestroyArray(&node);
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  n_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *form, const
  char_T *identifier, char_T y[2])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  j_emlrt_marshallIn(sp, emlrtAlias(form), &thisId, y);
  emlrtDestroyArray(&form);
}

static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[2])
{
  o_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[2]
{
  real_T (*ret)[2];
  int32_T iv6[2];
  int32_T i3;
  for (i3 = 0; i3 < 2; i3++) {
    iv6[i3] = 1 + i3;
  }

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, iv6);
  ret = (real_T (*)[2])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static void m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  int32_T iv7[1];
  boolean_T bv0[1] = { true };

  static const int32_T iv8[1] = { -1 };

  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 1U, iv8, &bv0[0],
    iv7);
  ret->size[0] = iv7[0];
  ret->allocatedSize = ret->size[0];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static void n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  int32_T iv9[2];
  int32_T i4;
  int32_T iv10[2];
  boolean_T bv1[2] = { true, false };

  for (i4 = 0; i4 < 2; i4++) {
    iv9[i4] = 3 * i4 - 1;
  }

  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 2U, iv9, &bv1[0],
    iv10);
  ret->size[0] = iv10[0];
  ret->size[1] = iv10[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[2])
{
  int32_T iv11[2];
  int32_T i5;
  for (i5 = 0; i5 < 2; i5++) {
    iv11[i5] = 1 + i5;
  }

  emlrtCheckBuiltInR2012b(sp, msgId, src, "char", false, 2U, iv11);
  emlrtImportCharArrayR2014b(sp, src, ret, 2);
  emlrtDestroyArray(&src);
}

void mlsLinearBasis2D_api(const mxArray * const prhs[5], const mxArray *plhs[1])
{
  emxArray_real_T *b_index;
  emxArray_real_T *node;
  emxArray_real_T *di;
  emxArray_real_T *phi;
  real_T (*pt)[2];
  char_T form[2];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real_T(&st, &b_index, 1, &c_emlrtRTEI, true);
  b_emxInit_real_T(&st, &node, 2, &c_emlrtRTEI, true);
  emxInit_real_T(&st, &di, 1, &c_emlrtRTEI, true);
  emxInit_real_T(&st, &phi, 1, &c_emlrtRTEI, true);

  /* Marshall function inputs */
  pt = c_emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "pt");
  e_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "index", b_index);
  g_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "node", node);
  e_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "di", di);
  i_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "form", form);

  /* Invoke the target function */
  mlsLinearBasis2D(&st, *pt, b_index, node, di, form, phi);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(phi);
  phi->canFreeData = false;
  emxFree_real_T(&phi);
  di->canFreeData = false;
  emxFree_real_T(&di);
  node->canFreeData = false;
  emxFree_real_T(&node);
  b_index->canFreeData = false;
  emxFree_real_T(&b_index);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_mlsLinearBasis2D_api.c) */
