function [psnr, ssim, fsim] = quality(imagery1, imagery2)
%==========================================================================
% Evaluates the quality assessment indices for two tensors.
%
% Syntax:
%   [psnr, ssim, fsim] = quality(imagery1, imagery2)
%
% Input:
%   imagery1 - the reference tensor
%   imagery2 - the target tensor

% NOTE: the tensor is a M*N*K array and DYNAMIC RANGE [0, 1]. 

% Output:
%   psnr - Peak Signal-to-Noise Ratio
%   ssim - Structure SIMilarity
%   fsim - Feature SIMilarity

% See also StructureSIM, FeatureSIM
%
% by Yi Peng
% update by Yu-Bang Zheng 11/19/2018
% update by Quan Yu 12/22/2021
%==========================================================================
imagery1=imagery1*255;imagery2=imagery2*255;
Nway = size(imagery1);
psnr = psnr_index(imagery1, imagery2);
fsim = FeatureSIM(imagery1, imagery2);
if ndims(imagery1) == 3 %images are colorful
    for i = 1:Nway(3)
        ssim(i) = ssim_index(imagery1(:, :, i), imagery2(:, :, i));
    end
    ssim = mean(ssim);
else
    ssim = ssim_index(imagery1, imagery2);
end


