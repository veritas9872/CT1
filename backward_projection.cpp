/*
 * Implementation of back-projection with no filtering.
 * Due to the difficulty of integrating FFT into a mex program,
 * I have decided to perform filtering in MATLAB.
 * Also, type checking is better done in MATLAB than in the mex program.
 * As such, type checking code has been removed.
 *
 * The thinking and terminology has been changed to match MATLAB's implementation.
 * The detector is on the line y*cos(\theta) = x*sin(\theta) + dist2det.
 * Indexing of the pixels is from left to right with the incoming ray as the top.
 * */


#include <mex.h>

#define _USE_MATH_DEFINES  // Necessary on older systems to access PI and other constants.
#include <cmath>

#include <iostream>
#include <vector>

void error_check(int nlhs, int nrhs, const mxArray **prhs);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    using std::cout; using std::endl;

    error_check(nlhs, nrhs, prhs);

    float *input_array;

    // Fetching input data.
#if MX_HAS_INTERLEAVED_COMPLEX
    input_array = mxGetSingles(prhs[0]);
#else
    input_array  = mxGetPr(prhs[0]);
#endif

    // Fetching input variables.
    const size_t num_img_rows = mxGetScalar(prhs[1]);
    const size_t num_img_cols = mxGetScalar(prhs[2]);
    const float img_pix_len_x = mxGetScalar(prhs[3]);
    const float img_pix_len_y = mxGetScalar(prhs[4]);
    const float det_pix_len = mxGetScalar(prhs[5]);
    const float projection_range = mxGetScalar(prhs[6]);

    const size_t num_det_pix = mxGetM(prhs[0]);
    const size_t num_views = mxGetN(prhs[0]);
    const size_t num_elements = mxGetNumberOfElements(prhs[0]);

    // Defining outputs.
    plhs[0] = mxCreateNumericMatrix(num_img_rows, num_img_cols, mxSINGLE_CLASS, mxREAL);

    float *output_array;
    output_array = mxGetSingles(plhs[0]);

    const std::vector<float> filtered_sinogram(input_array, input_array + num_elements);
    std::vector<float> reconstruction(num_img_rows * num_img_cols, 0);

    float theta;
    float x0, y0, x1, y1;
    float cos_th, sin_th, cos_th_sq, sin_th_sq, cos_th_sin_th;
    float det_center_x, det_center_y, det_start_x, det_start_y;
    float dist_from_start, pix_from_start, det_idx_delta;
    size_t det_idx;
    bool in_range;
    float linear_interp;
    const float radian = M_PIf32 / 180;  // (PI/180), used for degree to radian conversion.
    const float radian_delta = projection_range / num_views * radian;
    const float img_len_x = img_pix_len_x * num_img_cols;
    const float img_len_y = img_pix_len_y * num_img_rows;
    const float img_offset_x = static_cast<float>(num_img_cols - 1) * img_pix_len_x / 2;
    const float img_offset_y = static_cast<float>(num_img_rows - 1) * img_pix_len_y / 2;
    const float det_offset = static_cast<float>(num_det_pix - 1) * det_pix_len / 2;
    // Arbitrary assignment of distance to detector. Not important in parallel-beam geometry.
    const float dist2det = sqrtf(img_len_x * img_len_x + img_len_y * img_len_y + 3);
    const float dist2offset_sq = dist2det * dist2det + det_offset * det_offset;
    // scaling for xy axis is actually wrong.
    // This code only works if the pixel lengths for the x and y axes are the same.
    // TODO: Figure out why 1/2 is necessary in scaling!
    const float scale = 0.5f * M_PIf32 / num_views * sqrtf(img_pix_len_x * img_pix_len_y) / det_pix_len;

    std::vector<float> thetas(num_views, 0);
    for (int idx = 0; idx < num_views; idx++) {
        thetas.at(idx) = static_cast<float>(idx) * radian_delta;
    }
    std::vector<float> x0s(num_img_cols, 0);
    for (int idx = 0; idx < num_img_cols; idx++) {
        x0s.at(idx) = static_cast<float>(idx) * +img_pix_len_x - img_offset_x;
    }
    std::vector<float> y0s(num_img_rows, 0);
    for (int idx = 0; idx < num_img_rows; idx++) {
        y0s.at(idx) = static_cast<float>(idx) * -img_pix_len_y + img_offset_y;
    }

    for (int view = 0; view < num_views; view++) {
        theta = thetas.at(view);
        cos_th = cosf(theta);
        sin_th = sinf(theta);
        cos_th_sq = cos_th * cos_th;
        sin_th_sq = sin_th * sin_th;
        cos_th_sin_th = cos_th * sin_th;
        det_center_x = dist2det * -sin_th;  // cos(\theta + \pi/2) = -sin(\theta)
        det_center_y = dist2det * +cos_th;  // sin(\theta + \pi/2) = +cos(\theta)
        det_start_x = det_center_x - det_offset * cos_th;
        det_start_y = det_center_y -  det_offset * sin_th;


        for (int col = 0; col < num_img_cols; col++) {  // Column-major access for MATLAB data.
            // Not sure if it is better to pre-compute and store the x,y coordinates in a vector
            // or just re-calculate them again each time. It may depend on input size as well.

            for (int row = 0; row < num_img_rows; row++) {
                x0 = x0s.at(col); y0 = y0s.at(row);
                x1 = x0 * cos_th_sq + y0 * cos_th_sin_th - dist2det * sin_th;
                y1 = y0 * sin_th_sq + x0 * cos_th_sin_th + dist2det * cos_th;

                // Checks if the point is inside the detector range.
                // Excludes the edges of the detector.
                // Although a complicated interpolation might be more accurate,
                // a simpler solution would be to zero-pad the sinogram
                // in the detector (height) dimension by 1 pixel.
                in_range = (x1 * x1 + y1 * y1) < dist2offset_sq;

                if (in_range) {
                    dist_from_start = sqrtf((x1 - det_start_x) * (x1 - det_start_x)
                                            + (y1 - det_start_y) * (y1 - det_start_y));

                    pix_from_start = dist_from_start / det_pix_len;

                    // Converting float to size_t automatically floors the value too.
                    det_idx = floorf(pix_from_start);
                    det_idx_delta = pix_from_start - det_idx;

                    linear_interp = filtered_sinogram.at(det_idx + view * num_det_pix) * (1 - det_idx_delta)
                                    + filtered_sinogram.at(det_idx + 1 + view * num_det_pix) * det_idx_delta;

                    // Assign scaled value to reconstruction location.
                    reconstruction.at(row + col * num_img_rows) += linear_interp * scale;
                }
            }
        }
    }
    std::copy(reconstruction.begin(), reconstruction.end(), output_array);
}

// Error checking code.
void error_check(int nlhs, int nrhs, const mxArray **prhs){
    /* Check for proper number of arguments. */
    if(nrhs!=7) {
        mexErrMsgIdAndTxt("CT1:backward_projection:nrhs","Seven inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("CT1:backward_projection:nlhs","One output required.");
    }
    /* Make sure the first input argument is a single-type matrix. */
    if( !mxIsSingle(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("CT1:backward_projection:notSingle","Input must be a single-type floating point array.");
    }
    /* Check that the input is a matrix. */
    if(mxGetNumberOfDimensions(prhs[0])!=2) {
        mexErrMsgIdAndTxt("CT1:backward_projection:notMatrix","The first input must be a 2D matrix.");
    }
    /* Check that all other inputs are scalars. */
    for (size_t n=1; n < nrhs; n++) {
        if( !mxIsScalar(prhs[n]) || mxIsComplex(prhs[n]) ) {
            mexErrMsgIdAndTxt("CT1:backward_projection:notRealScalar",
                              "All inputs except the first must be (real-valued) scalars.");
        }
    }
}
