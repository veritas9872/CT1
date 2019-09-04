#include <mex.h>

#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <vector>


void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]);


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
    size_t num_det_pix = mxGetScalar(prhs[1]);  // Number of detector pixels.
    float det_pix_len = mxGetScalar(prhs[2]);  // Detector pixel length.

    size_t num_img_pix_x = mxGetScalar(prhs[3]);
    size_t num_img_pix_y = mxGetScalar(prhs[4]);

    float img_pix_len_x = mxGetScalar(prhs[5]);
    float img_pix_len_y = mxGetScalar(prhs[6]);

    float sampling_interval = mxGetScalar(prhs[7]);
    size_t num_views = mxGetScalar(prhs[8]);
    float projection_range = mxGetScalar(prhs[9]);  // Should be in degrees, not radians.

    if ((img_pix_len_x < sampling_interval) || (img_pix_len_y < sampling_interval)) {
        cout << "Warning! Sampling interval is shorter than pixel length!" << endl;
    }

    /* get dimensions of the input matrix */
    size_t num_cols = mxGetN(prhs[0]);
    size_t num_rows = mxGetM(prhs[0]);
    size_t num_elements = mxGetNumberOfElements(prhs[0]);

    cout << "Num cols: " << num_cols << ". Num rows: " << num_rows << ". Total: " << num_elements << '.' << endl;

    // Unavoidable memory copy. No way to transfer data from C array to std vector without memory copy.
    // If memory copy can be removed while vector is still used, please do so. Memory copy is very time consuming.
    std::vector<float> image(input_array, input_array + num_elements);

    /* create the output matrix */
    // Sinogram dimensions. M: rows = num_det_pix, N: cols = num_views
    plhs[0] = mxCreateNumericMatrix(num_det_pix, num_views, mxSINGLE_CLASS, mxREAL);

    /* get a pointer to the real data in the output matrix */
#if MX_HAS_INTERLEAVED_COMPLEX
    output_array = mxGetSingles(plhs[0]);
#else
    output_array = mxGetPr(plhs[0]);
#endif

    float img_len_x = num_img_pix_x * img_pix_len_x;
    float img_len_y = num_img_pix_y * img_pix_len_y;
    float det_len = num_det_pix * det_pix_len;

    if ((det_len < img_len_x) || (det_len < img_len_y)) {
        cout << "Warning! Detector is too short to cover the entire image!" << endl;
    }

    // Distance from center to source and device have been chosen arbitrarily.
    // They are not very important in parallel-beam geometry.
    float dist2src = sqrtf(img_len_x * img_len_x + img_len_y * img_len_y + 3);
    float dist2dst = dist2src;

    size_t num_acquisitions = floorf((dist2src + dist2dst) / sampling_interval);

    // Computational Routine
    float phi;
    float rad = (M_PIf32 / 180);  // (PI / 180) in float.
    float radian_delta = projection_range / num_views * rad;

    float img_range_x = (num_img_pix_x - 1) * img_pix_len_x / 2;
    float img_range_y = (num_img_pix_y - 1) * img_pix_len_y / 2;

    float det_offset = (num_det_pix - 1) * det_pix_len / 2;
    float det_offset_x, det_offset_y;

    float x_delta, y_delta;
    float x_center, y_center;
    float x_begin, y_begin;
    float x_pos, y_pos;
    float x_shift, y_shift;

    size_t view, det_pix_idx, acq;
    bool in_range;
    float ray_sum;

    float top_left, top_right, down_left, down_right;
    float dc, dr;
    size_t row, col;
    float bilinear;

    std::vector<float> sinogram(num_views * num_det_pix, 0);

    for (view = 0; view < num_views; view++) {
        // Rotating clockwise with the source at the top and detector at the bottom in the beginning.
        phi = view * radian_delta;  // Angle between detector and x-axis.

        // xy coordinates of the source center.
        x_center = dist2src * sinf(phi);
        y_center = dist2src * cosf(phi);

        // Distance from center of detector to fist/last pixel.
        det_offset_x = det_offset * cosf(phi);
        det_offset_y = det_offset * sinf(phi);

        // Distance between two adjacent pixels.
        x_shift = det_pix_len * cosf(phi);
        y_shift = det_pix_len * sinf(phi);

        // Movement of X-ray in each direction for each sampling interval.
        x_delta = -sinf(phi) * sampling_interval;
        y_delta = -cosf(phi) * sampling_interval;

        for (det_pix_idx = 0; det_pix_idx < num_det_pix; det_pix_idx++) {
            // The source has numbering from 0~n-1 from the left when the X-ray projection direction is up.
            // At the start, this means that 0 is on the right of the source and n-1 is on the left of the source.
            x_begin = x_center + (det_offset_x - det_pix_idx * x_shift);
            y_begin = y_center - (det_offset_y - det_pix_idx * y_shift);

            ray_sum = 0;

            for (acq = 0; acq < num_acquisitions; acq++) {
                x_pos = x_begin + acq * x_delta;
                y_pos = y_begin + acq * y_delta;

                // Ignoring cases where X-rays hit the edges. This will be fixed later.
                // Alternatively, the input data could be zero-padded to allow this code to work exactly.
                in_range = (-img_range_x <= x_pos) && (x_pos < img_range_x)
                           && (-img_range_y < y_pos) && (y_pos <= img_range_y);

                if (in_range) {
                    // Changing from xy coordinates to row/column coordinate system of the input image.
                    col = floorf(img_range_x + x_pos);
                    row = floorf(img_range_y - y_pos);
                    dc = img_range_x + x_pos - col;
                    dr = img_range_y - y_pos - row;

                    // Column major indexing of the input.
                    top_left = image.at(col * num_rows + row);
                    down_left = image.at(col * num_rows + (row + 1));
                    top_right = image.at((col + 1) * num_rows + row);
                    down_right = image.at((col + 1) * num_rows + (row + 1));

                    bilinear = top_left * (1 - dr) * (1 - dc)
                               + down_left * dr * (1 - dc)
                               + top_right * (1 - dr) * dc
                               + down_right * dr * dc;

                    ray_sum += bilinear;
                }
            }
            // "sinogram" should be in column-major order.
            sinogram.at(view * num_det_pix + det_pix_idx) += ray_sum;
        }
    }
    // Unavoidable memory copy again. Please remove this if it becomes possible at a later date.
    std::copy(sinogram.begin(), sinogram.end(), output_array);
}


// Error checking code.
void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]){
    /* check for proper number of arguments */
    if(nrhs!=10) {
        mexErrMsgIdAndTxt("CT1:forward_projection:nrhs","Ten inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("CT1:forward_projection:nlhs","One output required.");
    }
    /* make sure the first input argument is scalar */
    if( !mxIsSingle(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("CT1:forward_projection:notSingle","Input must be a single-type floating point array.");
    }

    /* check that number of rows in second input argument is 1 */
    if(mxGetNumberOfDimensions(prhs[0])!=2) {
        mexErrMsgIdAndTxt("CT1:forward_projection:notMatrix","The first input must be a 2D matrix.");
    }

    /* Checking size*/
    size_t num_cols = mxGetN(prhs[0]);
    size_t num_rows = mxGetM(prhs[0]);
    size_t num_elements = mxGetNumberOfElements(prhs[0]);
    if (num_elements != num_cols * num_rows){
        mexErrMsgIdAndTxt("CT1:forward_projection:wrongShape","Impossible shape for a 2D matrix.");
    }

    // Maybe add type checks for the others too.
    for (int n=1; n < nrhs; n++) {
        if(!mxIsScalar(prhs[n])) {
            mexErrMsgIdAndTxt("CT1:forward_projection:notScalar","All inputs except the first must be scalars.");
        }
    }
}