function [X,Y,V] = DTRTC_solver(data,known,Nway,coreNway,opts)
%% Data preprocessing and initialization
for u=1:2
    n1{u}=Nway{u}(1);n2{u}=Nway{u}(2);n3{u}=Nway{u}(3);
end
rank_min = opts.rank_min;
normM=norm(opts.Mtr(:));
maxit=opts.maxIter;
TC=opts.TC;
[known,id] =sort(known);
data =data(id);
gamma=opts.gamma;
t(1)=1;t(2)=1;
%% initialization X{u} and Y{u}
for u=1:2
    TC{u}=fft(TC{u},[],3);
    for i = 1:ceil((n3{u}+1)/2)
        if u<=2
            [Uk,Sigmak,Vk]=svds(TC{u}(:,:,i),coreNway{u}(i));
        else
            Uk=randn(n1{u},coreNway{u}(i));Sigmak=diag(randn(coreNway{u}(i),1));Vk=randn(n2{u},coreNway{u}(i));
        end
    X{u,i} = Uk*Sigmak;
    Y{u,i} = Vk';
    end
end
%% compute the initialization residual
for u=1:2
    C{u}=zeros(Nway{u});
    for n = 1:ceil((n3{u}+1)/2)
        C{u}(:,:,n)=X{u,n}*Y{u,n};
    end
    for n = ceil((n3{u}+1)/2)+1:n3{u}
        C{u}(:,:,n)=conj(C{u}(:,:,n3{u}-n+2));
    end
    TC2{u}=ifft(C{u},[],3);
end
TTC=Fold(Unfold(TC2{2},Nway{2},1),Nway{1},3);
TC=(TC2{1}+gamma*TTC)/(1+gamma);
TC(known) = data;TC=real(TC);
C{1}=fft(TC,[],3);C{2}=fft(Fold(Unfold(TC,Nway{1},3),Nway{2},1),[],3);
rho=0.95; mu=1.01; p=1/2;
for k = 1:maxit
%% update (X{u},Y{u})
  oldTC=TC;
  [gamma1,gamma2]=Lp(norm(TC2{1}(:)-TC(:),2)^2,norm(TTC(:)-TC(:),2)^2,p);
%   [gamma1,gamma2]=Geman(norm(TC2{1}(:)-TC(:),2)^2,norm(TTC(:)-TC(:),2)^2,p);
  if k<=0
      %--------- X---------------------
      for n = 1:ceil((n3{1}+1)/2)
          X{1,n}=C{1}(:,:,n)*Y{1,n}'*pinv(Y{1,n}*Y{1,n}');
          Xsq{1,n} = X{1,n}'*X{1,n};
          C{1}(:,:,n)=X{1,n}*Y{1,n};
      end
      for n = ceil((n3{1}+1)/2)+1:n3{1}
          C{1}(:,:,n)=conj(C{1}(:,:,n3{1}-n+2));
      end
      TC2{1}=ifft(C{1},[],3);TTC=Fold(Unfold(TC2{2},Nway{2},1),Nway{1},3);
      gamma1=p*(norm(TC2{1}(:)-TC(:),2))^(2*p-2); gamma2=p*(norm(TTC(:)-TC(:),2))^(2*p-2);
      TC=(gamma1*TC2{1}+gamma2*TTC)/(gamma1+gamma2);
      TC=real(TC);TC(known)=data;%TC=max(TC,0);TC=min(TC,1);
      TC=TC+(t(k)-1)/t(k+1)*(TC-oldTC);oldTC=TC;
      C{1}=fft(TC,[],3);
      for n = 1:ceil((n3{1}+1)/2)
          Y{1,n} = pinv(Xsq{1,n})*X{1,n}'*C{1}(:,:,n);
          C{1}(:,:,n)=X{1,n}*Y{1,n};
      end
      for n = ceil((n3{1}+1)/2)+1:n3{1}
          C{1}(:,:,n)=conj(C{1}(:,:,n3{1}-n+2));
      end
      TC2{1}=ifft(C{1},[],3);TTC=Fold(Unfold(TC2{2},Nway{2},1),Nway{1},3);
      gamma1=p*(norm(TC2{1}(:)-TC(:),2))^(2*p-2); gamma2=p*(norm(TTC(:)-TC(:),2))^(2*p-2);
      TC=(gamma1*TC2{1}+gamma2*TTC)/(gamma1+gamma2);
      TC=real(TC);TC(known)=data;%TC=max(TC,0);TC=min(TC,1);
      TC=TC+(t(k)-1)/t(k+1)*(TC-oldTC);oldTC=TC;
      %------------------Y-------------------------
      C{2}=fft(Fold(Unfold(TC,Nway{1},3),Nway{2},1),[],3);
      for n = 1:ceil((n3{2}+1)/2)
          X{2,n}=C{2}(:,:,n)*Y{2,n}'*pinv(Y{2,n}*Y{2,n}');
          Xsq{2,n} = X{2,n}'*X{2,n};
          C{2}(:,:,n)=X{2,n}*Y{2,n};
      end
      for n = ceil((n3{2}+1)/2)+1:n3{2}
          C{2}(:,:,n)=conj(C{2}(:,:,n3{2}-n+2));
      end
      TC2{2}=ifft(C{2},[],3); TTC=Fold(Unfold(TC2{2},Nway{2},1),Nway{1},3);
      gamma1=p*(norm(TC2{1}(:)-TC(:),2))^(2*p-2); gamma2=p*(norm(TTC(:)-TC(:),2))^(2*p-2);
      TC=(gamma1*TC2{1}+gamma2*TTC)/(gamma1+gamma2);
      TC=real(TC);TC(known)=data;%TC=max(TC,0);TC=min(TC,1);
      TC=TC+(t(k)-1)/t(k+1)*(TC-oldTC);oldTC=TC;
      C{2}=fft(Fold(Unfold(TC,Nway{1},3),Nway{2},1),[],3);
      for n = 1:ceil((n3{2}+1)/2)
          Y{2,n} = pinv(Xsq{2,n})*X{2,n}'*C{2}(:,:,n);
          C{2}(:,:,n)=X{2,n}*Y{2,n};
      end
  else
      C{2}=fft(Fold(Unfold(TC,Nway{1},3),Nway{2},1),[],3);
      for u=1:2
          for n = 1:ceil((n3{u}+1)/2)
              X{u,n}=C{u}(:,:,n)*Y{u,n}'*pinv(Y{u,n}*Y{u,n}');
              Xsq{u,n} = X{u,n}'*X{u,n};
              Y{u,n} = pinv(Xsq{u,n})*X{u,n}'*C{u}(:,:,n);
              C{u}(:,:,n)=X{u,n}*Y{u,n};
          end
      end
  end  
    %% adjust the rank of (X,Y)
for u=1:2
    if rho<1
        max_k=max(coreNway{u}(1:floor((n3{u}+1)/2)));
        sum_k=sum(coreNway{u}(1:floor((n3{u}+1)/2)));
        sigmas=zeros(max_k*(floor((n3{u}+1)/2)),1);
        for i=1:floor((n3{u}+1)/2)
            s = svd(Xsq{u,i});
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
            if(iidx>sum(rank_min{u}(length(floor((n3{u}+1)/2)))))  
                idx=floor((id(iidx+1:sum_k)-1)/max_k);
                for n=1:floor((n3{u}+1)/2)
                    num=length(find(idx==n-1));
                    if(num>0)
                        if coreNway{u}(n)-num>rank_min{u}(n)
                            coreNway{u}(n) = coreNway{u}(n)-num;
                        else
                            coreNway{u}(n) = rank_min{u}(n);
                        end
                        [Qx,Rx] = qr(X{u,n},0);
                        [Qy,Ry] = qr(Y{u,n}',0);
                        [U,S,V] = svd(Rx*Ry');
                        sigv = diag(S);
                        X{u,n} = Qx*U(:,1:coreNway{u}(n))*spdiags(sigv(1:coreNway{u}(n)),0,coreNway{u}(n),coreNway{u}(n));
                        Y{u,n} = (Qy*V(:,1:coreNway{u}(n)))';
                        C{u}(:,:,n)=X{u,n}*Y{u,n};
                    end
                end
            end
            rho=rho*mu;
        end
    end
end
%%
for u=1:2
    for n = ceil((n3{u}+1)/2)+1:n3{u}
        C{u}(:,:,n)=conj(C{u}(:,:,n3{u}-n+2));
    end
end
for u=1:2
    TC2{u}=ifft(C{u},[],3);
end
TTC=Fold(Unfold(TC2{2},Nway{2},1),Nway{1},3);
TC=(gamma1*TC2{1}+gamma2*TTC)/(gamma1+gamma2);
TC=real(TC);TC(known)=data;
if k<=2 TC=max(TC,0);TC=min(TC,1); end
%% judge whether converges 
res=norm(TC(:)-oldTC(:))/norm(oldTC(:));
if res<opts.tol || k==maxit  % 5e-3  1e-2
    TC(known) = data;
    TC=real(TC);
    V=TC;
    break;
end
res1(k)=norm(TC(:)-opts.Mtr(:))/normM;
TC=TC+(t(k)-1)/t(k+1)*(TC-oldTC);
t(k+2)=1/2+sqrt(1+4*(t(k+1))^2)/2;
%% update C{u}
C{1}=fft(TC,[],3);
end
TC=max(TC,0);TC=min(TC,1);V=TC;