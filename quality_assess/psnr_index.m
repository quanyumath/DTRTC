function psnr = psnr_index(x,y)
ndim =size(x); pr = prod(ndim,'all');
psnr=10*log10(pr*255^2/norm(x(:)-y(:))^2);
    




