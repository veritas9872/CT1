#include <mex.h>

#include <iostream>


// Error checking code.
void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]){
    /* check for proper number of arguments */
    if(nrhs!=10) {
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
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    using namespace std;  // Use std only inside the mexFunction, not globally.

    // Expected inputs. Arranged in order.
    float *input_array;

    size_t num_det_pix;
    float det_pix_len;

    size_t num_img_pix_x;
    size_t num_img_pix_y;

    float img_pix_len_x;
    float img_pix_len_y;

    float sampling_interval;
    size_t num_views;
    float projection_range;


    // Expected outputs.
    float *output_array;

    errorCheck(nlhs, nrhs, prhs);  // Checks for errors in inputs and outputs.

    /* create a pointer to the real data in the input matrix  */
#if MX_HAS_INTERLEAVED_COMPLEX
    input_array = mxGetSingles(prhs[0]);
#else
    input_array = mxGetPr(prhs[0]);
#endif

    /* get the value of the scalar input  */
    num_det_pix = mxGetScalar(prhs[1]);
    det_pix_len = mxGetScalar(prhs[2]);

    num_img_pix_x = mxGetScalar(prhs[3]);
    num_img_pix_y = mxGetScalar(prhs[4]);

    img_pix_len_x = mxGetScalar(prhs[5]);
    img_pix_len_y = mxGetScalar(prhs[6]);

    sampling_interval = mxGetScalar(prhs[7]);
    num_views = mxGetScalar(prhs[8]);
    projection_range = mxGetScalar(prhs[9]);


    /* get dimensions of the input matrix */
    size_t num_cols = mxGetN(prhs[0]);
    size_t num_rows = mxGetM(prhs[0]);

    /* create the output matrix */
    plhs[0] = mxCreateNumericMatrix(num_rows, num_cols, mxSINGLE_CLASS, mxREAL);

    /* get a pointer to the real data in the output matrix */
#if MX_HAS_INTERLEAVED_COMPLEX
    output_array = mxGetSingles(plhs[0]);
#else
    output_array = mxGetPr(plhs[0]);
#endif

    /* call the computational routine */
    // TODO: Work on the code from here. Get the other value. Type checks might also be added for the other values.
}