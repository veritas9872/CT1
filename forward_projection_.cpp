/*
 * Trying out the new MATLAB C++ mex API.
 * This requires MATLAB 2018a or above and C++11 or above with a 64bit CPU.
 * The new API is cleaner and has many new features.
 * It also integrates with C++ far better than the previous API.
 */

#include "mex.hpp"
#include "mexAdapter.hpp"

#define _USE_MATH_DEFINES  // Necessary on older systems to access PI and other constants.
#include <cmath>

#include <iostream>


using std::cout; using std::endl;

// MATLAB variable types
using matlab::mex::ArgumentList;
using matlab::data::Array;
using matlab::data::TypedArray;


class MexFunction : public matlab::mex::Function {
public:
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();  // Call MATLAB functions from C++
    matlab::data::ArrayFactory factory;  // Creates MATLAB arrays.

    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        check_arguments(outputs, inputs);
        // const keyword added to prevent accidental memory copying.
        const TypedArray<float> image = inputs[0];

        // Getting the values of scalar inputs. Indexing [0] is needed to get the raw data inside the Array class.
        const size_t num_det_pix = inputs[1][0];  // Number of detector pixels.
        const float det_pix_len = inputs[2][0];  // Detector pixel length.

        const float img_pix_len_x = inputs[3][0];
        const float img_pix_len_y = inputs[4][0];

        const float sampling_interval = inputs[5][0];
        const size_t num_views = inputs[6][0];
        const float projection_range = inputs[7][0];  // Should be in degrees, not radians.

        if ((img_pix_len_x < sampling_interval) || (img_pix_len_y < sampling_interval)) {
            cout << "Warning! Sampling interval is shorter than pixel length!" << endl;
        }

        // const std::vector<size_t> dims = inputs[0].getDimensions();  // Equivalent to version below.
        const matlab::data::ArrayDimensions dims = inputs[0].getDimensions();
        const size_t num_img_pix_rows = dims[0];
        const size_t num_img_pix_cols = dims[1];
        const size_t num_elements = inputs[0].getNumberOfElements();

        cout << "Num cols: " << num_img_pix_cols << ". Num rows: "
             << num_img_pix_rows << ". Total: " << num_elements << '.' << endl;

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
        // An offset is the distance from the center to the center of the first/last pixel.
        const float img_offset_x = static_cast<float>(num_img_pix_cols - 1) * img_pix_len_x / 2;
        const float img_offset_y = static_cast<float>(num_img_pix_rows - 1) * img_pix_len_y / 2;
        const float det_offset = static_cast<float>(num_det_pix - 1) * det_pix_len / 2;

        // Creating array that will be used as the output. A 2D matrix with (row: num_det_pix, col: num_views).
        TypedArray<float> sinogram = factory.createArray<float>({num_det_pix, num_views});

        // Computational Routine.
        // The for-loops here should parallelize very well as they are independent of one another.
        // Please take note in parallel implementations. Perhaps use multi-processing to parallelize the for-loops.
        for (size_t view = 0; view < num_views; view++) {
            // Rotating clockwise with the detector at the top and source at the bottom in the beginning.
            // Detector center begins at (0, dist2det) in xy coordinates and rotates clockwise.
            // This is different from the radon function in MATLAB.
            // Flip the output left-right to get equivalent outputs.
            phi = view * radian_delta;  // Angle between incoming X-ray and the y-axis.

            // xy coordinates of the detector center.
            det_center_x = dist2det * sinf(phi);
            det_center_y = dist2det * cosf(phi);

            // Movement of X-ray in each direction for each sampling interval.
            x_delta = sinf(phi) * sampling_interval;
            y_delta = cosf(phi) * sampling_interval;

            for (size_t det_pix_idx = 0; det_pix_idx < num_det_pix; det_pix_idx++) {
                // The source and detector have numbering from 0~n-1 from the left when the X-ray direction is up.
                // xy coordinates are in mm with the object center at the origin.
                // At the start, this means that 0 is on the right of the detector and n-1 is on the left.
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
                        // Changing from xy coordinates to row/column coordinate system of the input image.
                        // Each index value is set as the center of each pixel with that index with zero-indexing.
                        col = floorf((img_offset_x + x_pos) / img_pix_len_x);
                        row = floorf((img_offset_y - y_pos) / img_pix_len_y);
                        dc = (img_offset_x + x_pos) / img_pix_len_x - col;
                        dr = (img_offset_y - y_pos) / img_pix_len_y - row;

                        // Aligned and coalesced memory access is impossible, making this code memory inefficient.
                        // Indexing is done as (row, column) no matter what the underlying storage pattern is.
                        top_left = image[row][col];
                        down_left = image[row][col + 1];
                        top_right = image[row + 1][col];
                        down_right = image[row + 1][col + 1];

                        // Bilinear interpolation of nearby pixel values.
                        ray_sum += top_left * (1 - dc) * (1 - dr)
                                   + down_left * (1 - dc) * dr
                                   + top_right * dc * (1 - dr)
                                   + down_right * dc * dr;
                    }
                }
                // Indexing is done as (row, column) no matter what the underlying storage pattern is.
                // Column-major simply means that the last dimension of indexing is contiguous in memory.
                sinogram[det_pix_idx][view] = ray_sum;
            }
        }
        // Assign without memory copy through move semantics.
        outputs[0] = std::move(sinogram);
    }

    void check_arguments(ArgumentList outputs, ArgumentList inputs) {
        if (inputs.size() != 8) {
            matlabPtr->feval(u"error",0,std::vector<Array>(
                    { factory.createScalar("Eight inputs required") }));
        }
        if (outputs.size() != 1) {
            matlabPtr->feval(u"error",0,std::vector<Array>(
                    { factory.createScalar("One output required") }));
        }
        if (inputs[0].getType() != matlab::data::ArrayType::SINGLE ||
            inputs[0].getType() == matlab::data::ArrayType::COMPLEX_SINGLE ||
            inputs[0].getDimensions().size() != 2) {
            matlabPtr->feval(u"error",0,std::vector<Array>(
                    { factory.createScalar("First input must be a single-type matrix.") }));
        }
        for (size_t idx=1; idx < inputs.size(); idx++) {
            if (inputs[idx].getNumberOfElements() != 1) {
                matlabPtr->feval(u"error",0, std::vector<Array>(
                        { factory.createScalar("All inputs except the first must be scalars.") }));
            }
        }
        if (inputs[1].getDimensions().size() != 2) {
            matlabPtr->feval(u"error",0,std::vector<Array>(
                    { factory.createScalar("Input must be m-by-n dimension") }));
        }
    }
};