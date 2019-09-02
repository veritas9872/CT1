#include "mex.h"
#include <matrix.h>

#include <cassert>

#include <iostream>
#include <vector>


// Error checking code.
void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]){
    /* check for proper number of arguments */
    if(nrhs!=9) {
        mexErrMsgIdAndTxt("CT1:forward:nrhs","Ten inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("CT1:forward:nlhs","One output required.");
    }
    /* make sure the first input argument is scalar */
    if( !mxIsSingle(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("CT1:forward:notSingle","Input must be a single-type floating point array.");
    }

    /* check that number of rows in second input argument is 1 */
    if(mxGetNumberOfDimensions(prhs[0])!=2) {
        mexErrMsgIdAndTxt("CT1:forward:notMatrix","The first input must be a 2D matrix.");
    }

    // Maybe add type checks for the others too.
    for (int n=1; n < nrhs; n++) {
        if(!mxIsScalar(prhs[n])) {
            mexErrMsgIdAndTxt("CT1:forward:notScalar","All inputs except the first must be scalars.");
        }
    }
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    using std::cout; using std::endl;  // Use std only inside the mexFunction, not globally.

    // Expected inputs. Arranged in order.
    float *input_array;

    // Expected outputs.
    float *output_array;

    errorCheck(nlhs, nrhs, prhs);  // Checks for errors in inputs and outputs.

    /* create a pointer to the real data in the input matrix  */
    // Data comes in column-major as a C-array with 1 dimension.
#if MX_HAS_INTERLEAVED_COMPLEX
    input_array = mxGetSingles(prhs[0]);
#else
    input_array = mxGetPr(prhs[0]);
#endif

    /* get the value of the scalar input  */
    size_t num_det_pix = mxGetScalar(prhs[1]);
    float det_pix_len = mxGetScalar(prhs[2]);

    size_t num_img_pix_x = mxGetScalar(prhs[3]);
    size_t num_img_pix_y = mxGetScalar(prhs[4]);

    float img_pix_len_x = mxGetScalar(prhs[5]);
    float img_pix_len_y = mxGetScalar(prhs[6]);

    size_t num_views = mxGetScalar(prhs[7]);
    float projection_range = mxGetScalar(prhs[8]);

    /* get dimensions of the input matrix */
    size_t num_cols = mxGetN(prhs[0]);
    size_t num_rows = mxGetM(prhs[0]);
    size_t num_elements = mxGetNumberOfElements(prhs[0]);
    assert(num_elements == num_cols * num_rows);  // Checking size.

    // Unavoidable memory copy. No way to transfer data from C array to std vector without memory copy.
    // If memory copy can be removed while vector is still used, please do so. Memory copy is very time consuming.
    std::vector<float> input(input_array,input_array+num_elements);

    /* create the output matrix */
    plhs[0] = mxCreateNumericMatrix(num_rows, num_cols, mxSINGLE_CLASS, mxREAL);

    /* get a pointer to the real data in the output matrix */
#if MX_HAS_INTERLEAVED_COMPLEX
    output_array = mxGetSingles(plhs[0]);
#else
    output_array = mxGetPr(plhs[0]);
#endif

//    /* call the computational routine */
//    // TODO: Work on the code from here. Get the other value. Type checks might also be added for the other values.
    for (size_t idx = 0; idx < num_elements; idx++) {
        cout << input_array[idx] << endl;
    }
}