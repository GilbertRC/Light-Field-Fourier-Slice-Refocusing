% Copyright (c) 2020-  Richardson
% For research purpose only. Cannot be used for any other purpose without permission from the author(s).

% Inputs:
% -LF_FFT : Fourier transform of light field L(y,x,v,u,chan)(i.e. fftshift(fftn(ifftshift(LF))) ).
% -res_y : Resolution of y dimension.
% -res_x : Resolution of x dimension.
% -alpha : The refocusing depth ( disparity = 1 - 1/alpha ).
% -is_OS : Oversampling the slice or not.

% Outputs:
% -im_refocus : The refocused image.
% -img_fft : The Fourier transform of the refocused image.

function [im_refocus,im_fft] = refocus_fourierslice_rgb(LF_FFT, res_y, res_x, alpha, is_OS)
ny = res_y;
nx = res_x;
even_fft = 1-mod([ny nx],2);

[nyPad, nxPad, nvPad, nuPad, nChan] = size(LF_FFT);

% Calculate the frequency values for the spatial and angular domain
[kx,ky,kxMid,kyMid] = GenerateFrequencyGrid(nxPad,nyPad);
[ku,kv,kuMid,kvMid] = GenerateFrequencyGrid(nuPad,nvPad);

% Calculate the new frequency values for the angular domain through alpha
kv_alpha = (1 - alpha).*ky;
ku_alpha = (1 - alpha).*kx;

dKv = abs(kv(1,1) - kv(2,1)); % kv spacing
dKu = abs(ku(1,2) - ku(1,1)); % ku spacing

kv_alphaIdx = kvMid - (kv_alpha./dKv); % calculate the kv point based on pixel spacing for interpolation
% kv_alphaIdx( ( kv_alphaIdx - radius ) < 1 )= radius + 1; % Check for index values outside the valid range
% kv_alphaIdx( ( kv_alphaIdx ) > (length(kv) - radius )) = length(kv) - radius; % Check for index values outside the valid range
kv_alphaIdx( kv_alphaIdx < 1 )= 1;
kv_alphaIdx( kv_alphaIdx > length(kv) ) = length(kv);

ku_alphaIdx = (ku_alpha./dKu) + kuMid; % calculate the ku point based on pixel spacing for interpolation
% ku_alphaIdx( ( ku_alphaIdx - radius ) < 1 )= radius + 1; % Check for index values outside the valid range
% ku_alphaIdx( ( ku_alphaIdx ) > (length(ku) - radius )) = length(ku) - radius; % Check for index values outside the valid range
ku_alphaIdx( ku_alphaIdx < 1 )= 1;
ku_alphaIdx( ku_alphaIdx > length(ku) ) = length(ku);

[kxIdx, kyIdx] = meshgrid(1:nxPad,1:nyPad);

rIdx = ones(nyPad,nxPad,'double');
gIdx = 2*rIdx;
bIdx = 3*rIdx;

slice(:,:,1) = interpn(LF_FFT,kyIdx,kxIdx,kv_alphaIdx,ku_alphaIdx,rIdx,'linear');
slice(:,:,2) = interpn(LF_FFT,kyIdx,kxIdx,kv_alphaIdx,ku_alphaIdx,gIdx,'linear');
slice(:,:,3) = interpn(LF_FFT,kyIdx,kxIdx,kv_alphaIdx,ku_alphaIdx,bIdx,'linear');
slice(isnan(slice)) = 0;

%% Oversample the spectrum by 2 or not, perform ifft, then crop out central section
if is_OS == true
    kyIdx_OS = linspace(1,nyPad,2*nyPad);
    kxIdx_OS = linspace(1,nxPad,2*nxPad);
    [kxIdx_OS, kyIdx_OS] = meshgrid(kxIdx_OS, kyIdx_OS);
    kyMid_OS = ceil((2*nyPad+1)/2);
    kxMid_OS = ceil((2*nxPad+1)/2);
    slice_OS = interpn(kyIdx, kxIdx, slice, kyIdx_OS, kxIdx_OS);
    
    Temp = fftshift(ifftn(ifftshift(slice_OS)));
    im_refocus = real(Temp(kyMid_OS-floor(ny/2):kyMid_OS+floor(ny/2)-even_fft(1),kxMid_OS-floor(nx/2):kxMid_OS+floor(nx/2)-even_fft(2),:));
    im_fft = log(abs(slice_OS)+1);
    im_fft = im_fft(kyMid_OS-floor(ny/2):kyMid_OS+floor(ny/2)-even_fft(1),kxMid_OS-floor(nx/2):kxMid_OS+floor(nx/2)-even_fft(2),:);
else
    Temp = fftshift(ifftn(ifftshift(slice)));
    im_refocus = real(Temp(kyMid-floor(ny/2):kyMid+floor(ny/2)-even_fft(1),kxMid-floor(nx/2):kxMid+floor(nx/2)-even_fft(2),:));
    im_fft = log(abs(slice)+1);
    im_fft = im_fft(kyMid-floor(ny/2):kyMid+floor(ny/2)-even_fft(1),kxMid-floor(nx/2):kxMid+floor(nx/2)-even_fft(2),:);
end

end