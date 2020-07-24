% Copyright (c) 2020-  Richardson
% For research purpose only. Cannot be used for any other purpose without permission from the author(s).

clear all;

%% Parameters
res_y = 375;
res_x = 540;
res_v = 14;
res_u = res_v;
nChan = 3;
is_Padding = false;
is_OverSampling = false;

%% Load light field L(y,x,v,u,ch)
% if memory < 24G, please use gray image.
tic;fprintf('Load 4D Light Field...');
LF = zeros(res_y,res_x,res_v,res_u,nChan,'double');
for v=1:res_v
    for u=1:res_u
       img=im2double(imread(['Bikes_3_3/22/',num2str(v),'_',num2str(u),'.bmp']));
       LF(:,:,v,u,:)=img;
    end
end

t=toc;fprintf(['\b\b\b (done in ' num2str(t) 's)\n']);

%% Pre-processing (Pad+FFT)
tic;fprintf('Pre-processing for Fourier Slice Refocusing...');
% First, pad the LF before the FFT or not. 
% Good values are 5% in the spatial domain and 20 in the angular domain.
if is_Padding == true
    padVec = [round(0.05*res_y), round(0.05*res_x), 10, 10, 0];
    LF = padarray(LF, padVec); % zero pad the light field array
end

% Second, do FFT for the padded LF
LF_FFT = fftshift(fftn(ifftshift(LF)));

t=toc;fprintf(['\b\b\b (done in ' num2str(t) 's)\n']);

clear LF

%% Refocusing at alpha (Fourier Slice)
im_refocus = zeros(res_y,res_x,nChan,'double');
im_fft = zeros(res_y,res_x,nChan,'double');
i=1;
for alpha = 0.25:0.05:2
    tic;fprintf('Refocusing at alpha = %.2f',alpha);
    [im_refocus, im_fft] = refocus_fourierslice_rgb(LF_FFT, res_y, res_x, alpha, is_OverSampling);
    
    t=toc;fprintf([' (done in ' num2str(t) 's)\n']);
    
    % Saving the results of refocused image and its fourier transform
    scale_facter = max(max(max(im_refocus)));
    im_refocus=uint8(round(im_refocus/scale_facter*255));
    scale_facter = max(max(max(im_fft)));
    im_fft=uint8(round(im_fft/scale_facter*255));
    imwrite(im_fft,['fft_',num2str(i),'.bmp']);
    imwrite(im_refocus,[num2str(i),'.bmp']);
    i=i+1;
end
