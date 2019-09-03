%% 
% Code for using mex for CT.
% 
% Forward propagation with mex.

file_name = '/home/veritas/CLionProjects/CT1/forward_projection.cpp';

% mex file_name
mex -R2018a forward_projection.cpp

%%
input_array = single(magic(3));

num_det_pix=640;
det_pix_len=1;

num_img_pix_x=512;
num_img_pix_y=512;

img_pix_len_x=1;
img_pix_len_y=1;

sampling_interval=0.5;
num_views=360;
projection_range=360;
%%
my_sinogram = forward_projection(input_array, num_det_pix, det_pix_len,...
    num_img_pix_x, num_img_pix_y, img_pix_len_x, img_pix_len_y,...
    sampling_interval, num_views, projection_range);

%%






