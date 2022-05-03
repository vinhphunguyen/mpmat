/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_mlsLinearBasis2D_mex.c
 *
 * Code generation for function '_coder_mlsLinearBasis2D_mex'
 *
 */

/* Include files */
#include "mlsLinearBasis2D.h"
#include "_coder_mlsLinearBasis2D_mex.h"
#include "mlsLinearBasis2D_terminate.h"
#include "_coder_mlsLinearBasis2D_api.h"
#include "mlsLinearBasis2D_initialize.h"
#include "mlsLinearBasis2D_data.h"

/* Function Declarations */
static void mlsLinearBasis2D_mexFunction(int32_T nlhs, mxArray *plhs[1], int32_T
  nrhs, const mxArray *prhs[5]);

/* Function Definitions */
static void mlsLinearBasis2D_mexFunction(int32_T nlhs, mxArray *plhs[1], int32_T
  nrhs, const mxArray *prhs[5])
{
  int32_T n;
  const mxArray *inputs[5];
  const mxArray *outputs[1];
  int32_T b_nlhs;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 5) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 5, 4,
                        16, "mlsLinearBasis2D");
  }

  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 16,
                        "mlsLinearBasis2D");
  }

  /* Temporary copy for mex inputs. */
  for (n = 0; n < nrhs; n++) {
    inputs[n] = prhs[n];
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }

  /* Call the function. */
  mlsLinearBasis2D_api(inputs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  mlsLinearBasis2D_terminate();
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(mlsLinearBasis2D_atexit);

  /* Module initialization. */
  mlsLinearBasis2D_initialize();

  /* Dispatch the entry-point. */
  mlsLinearBasis2D_mexFunction(nlhs, plhs, nrhs, prhs);
}

/* End of code generation (_coder_mlsLinearBasis2D_mex.c) */
