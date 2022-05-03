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
#include "mlsLinearBasis2D_data.h"

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[2];
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_index,
  const char_T *identifier, real_T **y_data, int32_T y_size[1]);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T **y_data, int32_T y_size[1]);
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *node,
  const char_T *identifier))[2048];
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *pt, const
  char_T *identifier))[2];
static const mxArray *emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[1]);
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[2048];
static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *di,
  const char_T *identifier))[1024];
static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[1024];
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *form, const
  char_T *identifier, char_T y[12]);
static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[12]);
static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[2];
static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T **ret_data, int32_T ret_size[1]);
static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[2048];
static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[1024];
static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[12]);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[2]
{
  real_T (*y)[2];
  y = k_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_index,
  const char_T *identifier, real_T **y_data, int32_T y_size[1])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  d_emlrt_marshallIn(sp, emlrtAlias(b_index), &thisId, y_data, y_size);
  emlrtDestroyArray(&b_index);
}

static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T **y_data, int32_T y_size[1])
{
  l_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y_data, y_size);
  emlrtDestroyArray(&u);
}

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *node,
  const char_T *identifier))[2048]
{
  real_T (*y)[2048];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = f_emlrt_marshallIn(sp, emlrtAlias(node), &thisId);
  emlrtDestroyArray(&node);
  return y;
}
  static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *pt,
  const char_T *identifier))[2]
{
  real_T (*y)[2];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = b_emlrt_marshallIn(sp, emlrtAlias(pt), &thisId);
  emlrtDestroyArray(&pt);
  return y;
}

static const mxArray *emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[1])
{
  const mxArray *y;
  static const int32_T iv2[1] = { 0 };

  const mxArray *m2;
  y = NULL;
  m2 = emlrtCreateNumericArray(1, iv2, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m2, (void *)u_data);
  emlrtSetDimensions((mxArray *)m2, u_size, 1);
  emlrtAssign(&y, m2);
  return y;
}

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[2048]
{
  real_T (*y)[2048];
  y = m_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *di,
  const char_T *identifier))[1024]
{
  real_T (*y)[1024];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = h_emlrt_marshallIn(sp, emlrtAlias(di), &thisId);
  emlrtDestroyArray(&di);
  return y;
}

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[1024]
{
  real_T (*y)[1024];
  y = n_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *form,
  const char_T *identifier, char_T y[12])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  j_emlrt_marshallIn(sp, emlrtAlias(form), &thisId, y);
  emlrtDestroyArray(&form);
}

static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[12])
{
  o_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[2]
{
  real_T (*ret)[2];
  int32_T iv3[2];
  int32_T i1;
  for (i1 = 0; i1 < 2; i1++) {
    iv3[i1] = 1 + i1;
  }

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, iv3);
  ret = (real_T (*)[2])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T **ret_data, int32_T ret_size[1])
{
  int32_T iv4[1];
  boolean_T bv0[1] = { true };

  static const int32_T iv5[1] = { 52 };

  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 1U, iv5, &bv0[0],
    iv4);
  ret_size[0] = iv4[0];
  *ret_data = (real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
}

static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[2048]
{
  real_T (*ret)[2048];
  int32_T iv6[2];
  int32_T i2;
  for (i2 = 0; i2 < 2; i2++) {
    iv6[i2] = 1024 + -1022 * i2;
  }

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, iv6);
  ret = (real_T (*)[2048])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[1024]
{
  real_T (*ret)[1024];
  int32_T iv7[1];
  iv7[0] = 1024;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, iv7);
  ret = (real_T (*)[1024])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[12])
{
  int32_T iv8[2];
  int32_T i3;
  for (i3 = 0; i3 < 2; i3++) {
    iv8[i3] = 1 + 11 * i3;
  }

  emlrtCheckBuiltInR2012b(sp, msgId, src, "char", false, 2U, iv8);
  emlrtImportCharArrayR2014b(sp, src, ret, 12);
  emlrtDestroyArray(&src);
}

void mlsLinearBasis2D_api(const mxArray * const prhs[5], const mxArray *plhs[1])
{
  real_T (*phi_data)[52];
  real_T (*pt)[2];
  int32_T index_size[1];
  real_T (*index_data)[52];
  real_T (*node)[2048];
  real_T (*di)[1024];
  char_T form[12];
  int32_T phi_size[1];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  phi_data = (real_T (*)[52])mxMalloc(sizeof(real_T [52]));

  /* Marshall function inputs */
  pt = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "pt");
  c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "index", (real_T **)&index_data,
                     index_size);
  node = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "node");
  di = g_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "di");
  i_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "form", form);

  /* Invoke the target function */
  mlsLinearBasis2D(&st, *pt, *index_data, index_size, *node, *di, form,
                   *phi_data, phi_size);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*phi_data, phi_size);
}

/* End of code generation (_coder_mlsLinearBasis2D_api.c) */
