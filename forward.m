%% 
% Code for using mex for CT.
% 
% Forward propagation with mex.

file_name = '/home/veritas/CLionProjects/CT1/forward_projection.cpp';

% mex file_name
mex -R2018a forward_projection.cpp

%%
input_array = single(magic(3));

num_det_pix=3;
det_pix_len=1;  % In mm

img_pix_len_x=1;  % In mm
img_pix_len_y=1;  % In mm

sampling_interval=0.5;  % In mm
num_views=3;
projection_range=360;  % In degrees
%%
my_sinogram = forward_projection(input_array, num_det_pix, det_pix_len,...
    img_pix_len_x, img_pix_len_y, sampling_interval, num_views,...
    projection_range);

%%
imshow(my_sinogram)


