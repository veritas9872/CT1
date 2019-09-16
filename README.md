## CT Parallel-Beam Imaging

**This project contains the code for CT parallel beam imaging using MATLAB mex files.**

It contains code for forward projection and will soon include code for backward projection.

The code is designed for the MATLAB mex C and C++ APIs.

Files with the C++ API have an underscore at the end of their names.

The C API has better backward compatibility and documentation.

The C++ API integrates new C++11 features and requires MATLAB 2018a or newer.

All code assumes MATLAB version R2018a or higher and C++11 is used on a CPU with 64bit architecture .

A CMake file is included for easy use on IDEs.


## Known Issues.
1. The files cannot be compiled independently of the MATLAB mex function.
They can only be compiled in MATLAB using the mex function.
2. The C API does not allow direct memory allocation onto C++ std::vector objects. 
This requires additional memory copies and allocations, making the code inefficient.
3. The code has not been tested on older platforms.
4. On Ubuntu, the libstdc++.so.6. file may cause issues with compilation.
Change its name to libstdc++.so.6.old to solve compilation issues. 