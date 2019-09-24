%% 
% Code for using mex for CT.
% 
% Forward propagation with mex.
clear
mex -R2018a forward_projection.cpp
mex forward_projection_.cpp

%%
img_size = 256;
input_array = single(phantom(img_size));

num_det_pix=367;
det_pix_len=1;  % In mm

img_pix_len_x=1;  % In mm
img_pix_len_y=1;  % In mm

sampling_interval=1;  % In mm
num_views=270;
projection_range=360;  % In degrees

theta = (0:(num_views-1)) * (projection_range/num_views);

%%
sinogram = radon(input_array, theta);
figure(1)
imshow(sinogram, []); colorbar();

%%
my_sinogram = forward_projection(input_array, num_det_pix, det_pix_len,...
    img_pix_len_x, img_pix_len_y, sampling_interval, num_views,...
    projection_range);

%%
figure(2)
imshow(my_sinogram, []); colorbar();

sino_delta = sinogram - my_sinogram;
figure(3)
imshow(sino_delta, []); colorbar();

my_recon = iradon(my_sinogram, theta);
figure(4)
imshow(my_recon, []); colorbar();

recon = iradon(sinogram, theta);
figure(5)
imshow(recon, []); colorbar();

img_delta = my_recon - recon;
figure(6)
imshow(img_delta, []); colorbar();

%%
% Compare execution time with MATLAB function.

% orig = @() radon(input_array, theta);
% mine = @() forward_projection(input_array, num_det_pix, det_pix_len,...
%     img_pix_len_x, img_pix_len_y, sampling_interval, num_views,...
%     projection_range);
% cpp_mex = @() forward_projection_(input_array, num_det_pix, det_pix_len,...
%     img_pix_len_x, img_pix_len_y, sampling_interval, num_views,...
%     projection_range);
% 
% orig_time = timeit(orig);
% mine_time = timeit(mine);
% cpp_mex_time = timeit(cpp_mex);
% 
% disp("Time taken by MATLAB native function.")
% disp(orig_time)
% disp("Time taken by custom mex function.")
% disp(mine_time)
% disp("Time taken by the C++ mex API version.")
% disp(cpp_mex_time)
% 
% % A value lesser than 1 indicates slowdown
% disp("Relative Speedup for C mex")
% disp(orig_time/mine_time)
% disp("Relative Speedup for C++ mex")
% disp(orig_time/cpp_mex_time)
