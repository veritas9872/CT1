#include <mex.h>
#include <iostream>


/* The computational routine */


// Error checking code.
void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]){
    double dim;

    /* Check for proper number of arguments. */
    if ( nrhs != 1 || mxIsNumeric(prhs[0]) == 0 || mxGetNumberOfElements(prhs[0]) != 1 ) {
        mexErrMsgIdAndTxt("MATLAB:arraySize:rhs",
                          "This function requires one scalar numeric input.");
    }

    dim = mxGetScalar(prhs[0]);

    if ( dim < 0 ) {
        mexErrMsgIdAndTxt("MATLAB:arraySize:dimensionNegative",
                          "The input dimension must be positive.");
    }

#if !defined(MX_COMPAT_32)
    /* Make sure that it is safe to cast dim to mwSize when using largeArrayDims.*/
    if ( dim > MWSIZE_MAX ) {
        mexErrMsgIdAndTxt("MATLAB:arraySize:dimensionTooLarge",
                          "The input dimension, %.0f, is larger than the maximum value of mwSize, %u, when built with largeArrayDims.", dim, MWSIZE_MAX);
    }
#endif

    if ( nlhs > 1 ) {
        mexErrMsgIdAndTxt("MATLAB:arraySize:lhs","Too many output arguments.");
    }

}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // Expected inputs:
    size_t num_det_pix;
    float det_pix_len;

    size_t img_pix_num_x;
    size_t img_pix_num_y;

    float img_pix_len_x;
    float img_pix_len_y;

    float sampling_interval;
    size_t num_views;
    float projection_range;




    double multiplier;              /* input scalar */
    double *inMatrix;               /* 1xN input matrix */
    size_t ncols;                   /* size of matrix */
    double *outMatrix;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    /* make sure the first input argument is scalar */
    if( !mxIsDouble(prhs[0]) ||
        mxIsComplex(prhs[0]) ||
        mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input multiplier must be a scalar.");
    }

    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[1]) ||
        mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }

    /* check that number of rows in second input argument is 1 */
    if(mxGetM(prhs[1])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    }

    /* get the value of the scalar input  */
    multiplier = mxGetScalar(prhs[0]);

    /* create a pointer to the real data in the input matrix  */
#if MX_HAS_INTERLEAVED_COMPLEX
    inMatrix = mxGetDoubles(prhs[1]);
#else
    inMatrix = mxGetPr(prhs[1]);
#endif

    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[1]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)ncols,mxREAL);

    /* get a pointer to the real data in the output matrix */
#if MX_HAS_INTERLEAVED_COMPLEX
    outMatrix = mxGetDoubles(plhs[0]);
#else
    outMatrix = mxGetPr(plhs[0]);
#endif

    /* call the computational routine */
}