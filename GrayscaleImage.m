clc 
clear
addpath(genpath('TCTF-M'));
addpath(genpath('Data\Grayscale Image'));
addpath(genpath('quality_assess'));
%% 导入图片
indimgs = [1:30];
i=1;id = indimgs(i);
pic_name = [ 'Grayscale Image/',num2str(id),'.tiff'];
I = double(imread(pic_name));
XX = I/255;% imshow(XX)
%XX=rgb2gray(XX);
%% 样本率
sr=0.7;
Omega = find(rand(numel(XX),1)<sr);
%% 将矩阵重构为张量
k=16;  %  16  15
clear X
for i=1:k
    l=size(XX,2)/k;
    X(:,:,i)=XX(:,l*(i-1)+1:l*i);
end
%% TCTF-M
tic
TC=TCTFM(X,Omega,2);   
time_TCTFM=toc; 
for i=1:k
    TC_M(:,l*(i-1)+1:l*i)=TC(:,:,i);
end
[psnr_TCTFM,ssim_TCTFM,fsim_TCTFM]=quality(XX,TC_M);
%imshow(TC_M)   


