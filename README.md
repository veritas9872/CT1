## CT Parallel-Beam Imaging

**This project contains the code for CT parallel beam imaging using MATLAB mex files.**

It contains code for forward projection and will soon include code for backward projection.

The code is designed for the MATLAB mex C API, 
which has better compatibility and documentation than the newer C++ API.

It assumes MATLAB version R2018a or higher and C++11 is used.

A CMake file is included for easy use on IDEs.


## Known Issues.
1. The files cannot be compiled independently of the MATLAB mex function.
They can only be compiled in MATLAB.
2. The C API does not allow direct memory allocation onto C++ std::vector objects. 
This requires additional memory copies and allocations, making the code inefficient.