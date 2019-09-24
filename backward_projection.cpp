#include <mex.h>

#define _USE_MATH_DEFINES  // Necessary on older systems to access PI and other constants.
#include <cmath>

#include <iostream>
#include <vector>

#include <kfr>

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    using std::cout; using std::endl;

    const float *sinogram = mxGetSingles(prhs[0]);
    const size_t num_img_rows = mxGetScalar(prhs[1]);
    const size_t num_img_cols = mxGetScalar(prhs[2]);
    const float img_pix_len_x = mxGetScalar(prhs[3]);
    const float img_pix_len_y = mxGetScalar(prhs[4]);
    const float det_pix_len = mxGetScalar(prhs[5]);
    const float sampling_interval = mxGetScalar(prhs[6]);
    const float projection_range = mxGetScalar(prhs[7]);

    const size_t num_det_pix = mxGetM(prhs[0]);
    const size_t num_views = mxGetN(prhs[0]);



}
