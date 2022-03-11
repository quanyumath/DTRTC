function [X,Y,TC,TotalRes,Time] = TCTFM_time_solver(data,known,Nway,coreNway,opts,kk)
%% Data preprocessing and initialization
n1=Nway(1);n2=Nway(2);n3=Nway(3);maxit=opts.maxIter;tol=opts.tol;
[known,id] = sort(known);
data = data(id);
normM=norm(opts.Mtr(:));
rank_min = opts.rank_min;
%% initialization X and Y
TC=randn(Nway);
TC(known)=data;
TC=fft(TC,[],3);
for i = 1:floor((n3+1)/2)
    Uk=randn(n1,coreNway(i));Sigmak=diag(randn(coreNway(i),1));Vk=randn(n2,coreNway(i));
    X{i} = Uk*Sigmak;
    Y{i} = Vk';
   oldX{i}=X{i};oldY{i}=Y{i};
end

%% compute the initialization residual
TotalRes=[];Time=[];
C=zeros(n1,n2,n3);
for n = 1:floor((n3+1)/2)
    C(:,:,n)=X{n}*Y{n};
end
for n =floor((n3+1)/2)+1:n3
    C(:,:,n)=conj(C(:,:,n3-n+2));
end
TC=ifft(C,[],3);TC(known) = data;TC=real(TC);
C=fft(TC,[],3);

rho=0.95;
mu=1.01;
t(1)=1;t(2)=1;
%res=50;

for k = 1:maxit
    tic
     oldTC=TC;
    %% update (X,Y)
    if k<=kk
    for n = 1:floor((n3+1)/2)
        X{n}=C(:,:,n)*Y{n}'*pinv(Y{n}*Y{n}');X{n}=X{n}+(t(k)-1)/t(k+1)*(X{n}-oldX{n});oldX{n}=X{n};
        Xsq{n} = X{n}'*X{n};
        C(:,:,n)=X{n}*Y{n}; 
    end
    for n =floor((n3+1)/2)+1:n3
        C(:,:,n)=conj(C(:,:,n3-n+2));
    end
    TC=ifft(C,[],3);
    TC(known) = data;TC = max(TC,0);TC = min(TC,1);
    TC=TC+(t(k)-1)/t(k+1)*(TC-oldTC);oldTC=TC;
    C=fft(TC,[],3);
    for n = 1:floor((n3+1)/2)
        Y{n} =pinv(Xsq{n})*X{n}'*C(:,:,n);Y{n}=Y{n}+(t(k)-1)/t(k+1)*(Y{n}-oldY{n});oldY{n}=Y{n};
        C(:,:,n)=X{n}*Y{n};
    end
    else
    for n = 1:floor((n3+1)/2)
        X{n}=C(:,:,n)*Y{n}'*pinv(Y{n}*Y{n}');X{n}=X{n}+(t(k)-1)/t(k+1)*(X{n}-oldX{n});
        Xsq{n} = X{n}'*X{n};
        Y{n} =pinv(Xsq{n})*X{n}'*C(:,:,n);Y{n}=Y{n}+(t(k)-1)/t(k+1)*(Y{n}-oldY{n});
        C(:,:,n)=X{n}*Y{n};
        oldX{n}=X{n};oldY{n}=Y{n};
    end
    end
if 0
    %% adjust the rank of (X,Y)
    if rho<1
        max_k=max(coreNway(1:floor((n3+1)/2)));
        sum_k=sum(coreNway(1:floor((n3+1)/2)));
        sigmas=zeros(max_k*(floor((n3+1)/2)),1);
        for i=1:floor((n3+1)/2)
            s = svd(Xsq{i});
            sigmas((i-1)*max_k+1:(i-1)*max_k+length(s))=s;
        end
        [dR,id]=sort(sigmas,'descend');
        drops = dR(1:sum_k-1)./dR(2:sum_k);
        [dmx,imx] = max(drops);
        rel_drp = (sum_k-1)*dmx/(sum(drops)-dmx);
        if rel_drp>10
            thold=rho*sum(dR);
            iidx=0;ss=0;
            len=length(dR);
            for i=1:len
                ss=ss+dR(i);
                if(ss>thold)
                    iidx=i;
                    break;
                end
            end
            if(iidx>sum(rank_min(n)))
                idx=floor((id(iidx+1:sum_k)-1)/max_k);
                for n=1:floor((n3+1)/2)
                    num=length(find(idx==n-1));
                    if(num>0)
                        if coreNway(n)-num>rank_min(n)
                            coreNway(n) = coreNway(n)-num;
                        else
                            coreNway(n) = rank_min(n);
                        end
                        [Qx,Rx] = qr(X{n},0);
                        [Qy,Ry] = qr(Y{n}',0);
                        [U,S,V] = svd(Rx*Ry');
                        sigv = diag(S);
                        X{n} = Qx*U(:,1:coreNway(n))*spdiags(sigv(1:coreNway(n)),0,coreNway(n),coreNway(n));
                        Y{n} = (Qy*V(:,1:coreNway(n)))';
                        C(:,:,n)=X{n}*Y{n};
                    end
                end
            end
            rho=rho*mu;
        end
    end
end
    %% judge whether converges
    for n =floor((n3+1)/2)+1:n3
        C(:,:,n)=conj(C(:,:,n3-n+2));
    end
    TC=ifft(C,[],3);
%   res=norm(TC(known)-data)/normdata;
    TC=real(TC);
    TC(known) = data;TC = max(TC,0);TC = min(TC,1);
    res=norm(TC(:)-oldTC(:))/norm(oldTC(:));
    if res<tol 
        TC(known) = data;
        TC=real(TC);
        break;
    end
    if(k==maxit)
        break;
    end
     TC=TC+(t(k)-1)/t(k+1)*(TC-oldTC);
     t(k+2)=1/2+sqrt(1+4*(t(k+1))^2)/2;
    %% update C
    C=fft(TC,[],3);
   %% 统计每步指标
    % TotalRes(k)=norm(TC(:)-opts.Mtr(:))/normM;  % Numerical
    Time(k)=toc;sTime(k)=sum(Time(:));
    for i=1:15
        l=2250/15;
        TC_M(:,l*(i-1)+1:l*i)=TC(:,:,i);
        XX(:,l*(i-1)+1:l*i)=opts.Mtr(:,:,i);
    end    
    [psnr_TCTFM(k),ssim_TCTFM(k),fsim_TCTFM(k)]=quality(XX,TC_M);   %%  GrayscaleImage
end
    TotalRes=[psnr_TCTFM;ssim_TCTFM;fsim_TCTFM;sTime];
end