clc 
clear
addpath(genpath('DTRTC'));
addpath(genpath('Data\Color Image'));
addpath(genpath('quality_assess'));
%% 导入图片
indimgs = [1:10];
i=4;id = indimgs(i);
pic_name = [ 'Color Image/',num2str(id),'.tiff']; %test_image train
I = double(imread(pic_name));
X = I/255;% imshow(X)
%% 样本率
sr=0.4;
Omega = find(rand(numel(X),1)<sr);
G=zeros(size(X));G(Omega)=X(Omega);%imshow(G)
%% DTRTC
Nway=size(X);
Nway2=[Nway(3),Nway(2)*Nway(1)/50,50]; %重构张量的大小
tic
DT=DTRTC_Color(X,Omega,Nway2);   
time_DTRTC=toc;
[psnr_DTRTC,ssim_DTRTC,fsim_DTRTC]=quality(X,DT);
%imshow(DT)

