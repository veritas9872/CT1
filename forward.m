%% 
% Code for using mex for CT.
% 
% Forward propagation with mex.
clear
mex -R2018a forward_projection.cpp

%%
img_size = 256;
input_array = single(phantom(img_size));

num_det_pix=367;
det_pix_len=1;  % In mm

img_pix_len_x=1;  % In mm
img_pix_len_y=1;  % In mm

sampling_interval=1;  % In mm
num_views=180;
projection_range=180;  % In degrees

theta = 0:(projection_range/num_views):projection_range-1;
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

sino_delta = fliplr(sinogram) - my_sinogram;
figure(3)
imshow(sino_delta, []); colorbar();

my_recon = iradon(fliplr(my_sinogram), theta);
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
orig = @() radon(input_array, theta);
mine = @() forward_projection(input_array, num_det_pix, det_pix_len,...
    img_pix_len_x, img_pix_len_y, sampling_interval, num_views,...
    projection_range);

orig_time = timeit(orig);
mine_time = timeit(mine);

disp("Time taken by MATLAB native function.")
disp(orig_time)
disp("Time taken by custom mex function.")
disp(mine_time)
disp("Relative Speedup (A value lesser than 1 indicates slowdown.)")
disp(orig_time/mine_time)
