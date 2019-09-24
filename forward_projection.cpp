/*
 * Naive implementation forward projection with C++.
 * No attenuation factors or other details have been considered.
 * The code is not very efficient due to the use of std vectors.
 * This is because the MATLAB mex C API does not allow the direct allocation of MATLAB data to std vectors.
 * Instead, data must pass through C arrays.
 * This forces this code to allocate and copy to and from the C arrays to std vectors.
 * This results in two extra memory allocations and one extra zero-initialization.
 * Moreover, the data access pattern is not aligned and coalesced, making the program memory-bound, not compute bound.
 */

#include "mex.h"

#define _USE_MATH_DEFINES  // Necessary on older systems to access PI and other constants.
#include <cmath>

#include <iostream>
#include <vector>

void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]);

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    using std::cout; using std::endl;

    errorCheck(nlhs, nrhs, prhs);  // Checks for errors in inputs and outputs.

    float *input_array;
    float *output_array;

    /* create a pointer to the real data in the input matrix  */
    // Data comes in column-major as a C-array with 1 dimension.
#if MX_HAS_INTERLEAVED_COMPLEX
    input_array = mxGetSingles(prhs[0]);
#else
    input_array = mxGetPr(prhs[0]);
#endif

    /* get the value of the scalar input  */
    const size_t num_det_pix = mxGetScalar(prhs[1]);  // Number of detector pixels.
    const float det_pix_len = mxGetScalar(prhs[2]);  // Detector pixel length.

    const float img_pix_len_x = mxGetScalar(prhs[3]);
    const float img_pix_len_y = mxGetScalar(prhs[4]);

    const float sampling_interval = mxGetScalar(prhs[5]);
    const size_t num_views = mxGetScalar(prhs[6]);
    const float projection_range = mxGetScalar(prhs[7]);  // Should be in degrees, not radians.

    if ((img_pix_len_x < sampling_interval) || (img_pix_len_y < sampling_interval)) {
        cout << "Warning! Sampling interval is shorter than pixel length!" << endl;
    }

    /* get dimensions of the input matrix */
    const size_t num_img_pix_rows = mxGetM(prhs[0]);
    const size_t num_img_pix_cols = mxGetN(prhs[0]);
    const size_t num_elements = mxGetNumberOfElements(prhs[0]);

//    cout << "Num cols: " << num_img_pix_cols << ". Num rows: "
//         << num_img_pix_rows << ". Total: " << num_elements << '.' << endl;
//    assert(num_img_pix_rows * num_img_pix_cols == num_elements);

    // Unavoidable memory copy. No way to transfer data from C array to std vector without memory copy.
    // If memory copy can be removed while vector is still used, please do so. Memory copy is very time consuming.
    const std::vector<float> image(input_array, input_array + num_elements);

    /* create the output matrix */
    // Sinogram dimensions. M: rows = num_det_pix, N: cols = num_views
    plhs[0] = mxCreateNumericMatrix(num_det_pix, num_views, mxSINGLE_CLASS, mxREAL);

    /* get a pointer to the real data in the output matrix */
#if MX_HAS_INTERLEAVED_COMPLEX
    output_array = mxGetSingles(plhs[0]);
#else
    output_array = mxGetPr(plhs[0]);
#endif

    const float img_len_x = num_img_pix_cols * img_pix_len_x;
    const float img_len_y = num_img_pix_rows * img_pix_len_y;
    const float det_len = num_det_pix * det_pix_len;

    if ((det_len * det_len < img_len_x * img_len_y) || (det_len < img_len_x) || (det_len < img_len_y)) {
        cout << "Warning! Detector is too short to cover the entire image!" << endl;
    }
    if ((num_det_pix < num_img_pix_cols) || (num_det_pix < num_img_pix_rows)) {
        cout << "Warning! Detector pixel number is insufficient for image size!" << endl;
    }

    // Distance from center to source and device have been chosen arbitrarily.
    // They are not very important in parallel-beam geometry.
    const float dist2src = sqrtf(img_len_x * img_len_x + img_len_y * img_len_y + 3);  // Distance to source
    const float dist2det = dist2src;  // Distance to detector

    const size_t num_acquisitions = floorf((dist2src + dist2det) / sampling_interval);

    float phi;
    float det_center_x, det_center_y;
    float x_delta, y_delta;
    float x_end, y_end;
    float x_pos, y_pos;
    float dc, dr;
    size_t row, col;

    bool in_range;
    float ray_sum;
    float top_left, top_right, down_left, down_right;

    const float EPS = 1E-4;  // Necessary for numerical stability in edge cases of range.
    const float rad = (M_PIf32 / 180);  // (PI / 180) in float. Used for degree to radian conversion.
    const float radian_delta = projection_range / num_views * rad;
    // An offset is the distance from the center to the center of the last pixel.
    const float img_offset_x = static_cast<float>(num_img_pix_cols - 1) * img_pix_len_x / 2;
    const float img_offset_y = static_cast<float>(num_img_pix_rows - 1) * img_pix_len_y / 2;
    const float det_offset = static_cast<float>(num_det_pix - 1) * det_pix_len / 2;

    // Needed to replace division with multiplication.
    const float img_pix_scale_x = 1 / img_pix_len_x;
    const float img_pix_scale_y = 1 / img_pix_len_y;

    // Initializing sinogram with 0s. Unavoidable memory allocation by 0. Empty vector is not possible.
    std::vector<float> sinogram(num_det_pix * num_views, 0);

    // Computational Routine.
    // The for-loops here should parallelize very well as they are independent of one another.
    // Please take note of this in CUDA or OpenACC implementations.
    // Perhaps use multi-processing to parallelize the for-loops.
    for (size_t view = 0; view < num_views; view++) {
        // Changed code to make the behavior the same as MATLAB's radon function.
        // The detector starts at the bottom and rotates clockwise.
        // However, the code is now very confusing.
        phi = M_PIf32 - view * radian_delta;  // Angle between incoming X-ray and the y-axis.

        // xy coordinates of the detector center.
        det_center_x = dist2det * sinf(phi);
        det_center_y = dist2det * cosf(phi);

        // Movement of X-ray in each direction for each sampling interval.
        x_delta = sinf(phi) * sampling_interval;
        y_delta = cosf(phi) * sampling_interval;

        for (size_t det_pix_idx = 0; det_pix_idx < num_det_pix; det_pix_idx++) {
            // The source and detector have numbering from 0~n-1 from the left when the X-ray direction is up.
            // xy coordinates are in mm with the object center at the origin.
            // At the start, this means that 0 is on the right of the detector and n-1 is on the left of the detector.
            x_end = det_center_x + (det_offset - det_pix_idx * det_pix_len) * cosf(-phi);
            y_end = det_center_y + (det_offset - det_pix_idx * det_pix_len) * sinf(-phi);

            ray_sum = 0;  // Resetting the sum to 0 for each pixel.
            for (size_t acq = 0; acq < num_acquisitions; acq++) {
                // Going back from the detector to the source. Easier to code this way.
                x_pos = x_end - acq * x_delta;
                y_pos = y_end - acq * y_delta;

                // Ignoring cases where X-rays hit the edges. This will be fixed later.
                // Alternatively, the input data could be zero-padded to allow this code to work precisely.
                // EPS necessary for numerical stability.
                in_range = (-img_offset_x <= x_pos) && (x_pos < img_offset_x - EPS)
                           && (-img_offset_y + EPS < y_pos) && (y_pos <= img_offset_y);

                if (in_range) {
                    // Changing from xy coordinates to row/column coordinate system of the input image in column-major.
                    // Each index value is set as the center of each pixel with that index with zero-indexing.
                    col = floorf((img_offset_x + x_pos) * img_pix_scale_x);
                    row = floorf((img_offset_y - y_pos) * img_pix_scale_y);
                    dc = (img_offset_x + x_pos) * img_pix_scale_x - col;
                    dr = (img_offset_y - y_pos) * img_pix_scale_y - row;

//                    assert((0 <= dc) && (dc < 1) && (0 <= dr) && (dr < 1));

                    // Column major indexing of the input image.
                    // Aligned and coalesced memory access is impossible, making this code memory inefficient.
                    top_left = image.at(col * num_img_pix_rows + row);
                    down_left = image.at(col * num_img_pix_rows + (row + 1));
                    top_right = image.at((col + 1) * num_img_pix_rows + row);
                    down_right = image.at((col + 1) * num_img_pix_rows + (row + 1));

                    // Bilinear interpolation of nearby pixel values.
                    ray_sum += top_left * (1 - dc) * (1 - dr)
                            + down_left * (1 - dc) * dr
                            + top_right * dc * (1 - dr)
                            + down_right * dc * dr;
                }
            }
            // "sinogram" should be in column-major order.
            sinogram.at(view * num_det_pix + det_pix_idx) = ray_sum;
        }
    }
    // Unavoidable memory copy again. Please remove this if it becomes possible at a later date.
    std::copy(sinogram.begin(), sinogram.end(), output_array);
}


// Error checking code.
void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]){
    /* Check for proper number of arguments. */
    if(nrhs!=8) {
        mexErrMsgIdAndTxt("CT1:forward_projection:nrhs","Eight inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("CT1:forward_projection:nlhs","One output required.");
    }
    /* Make sure the first input argument is a single-type matrix. */
    if( !mxIsSingle(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("CT1:forward_projection:notSingle","Input must be a single-type floating point array.");
    }
    /* Check that the input is a matrix. */
    if(mxGetNumberOfDimensions(prhs[0])!=2) {
        mexErrMsgIdAndTxt("CT1:forward_projection:notMatrix","The first input must be a 2D matrix.");
    }
    /* Check that all other inputs are scalars. */
    for (size_t n=1; n < nrhs; n++) {
        if( !mxIsScalar(prhs[n]) || mxIsComplex(prhs[n]) ) {
            mexErrMsgIdAndTxt("CT1:forward_projection:notRealScalar",
                    "All inputs except the first must be (real-valued) scalars.");
        }
    }
}