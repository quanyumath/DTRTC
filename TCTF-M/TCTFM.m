function [TC,TotalRes,Time] = TCTFM( X,Omega,kk)
% kk:第一阶段迭代的步数

%% produce data
data=X(Omega);
known=Omega;
[n1,n2,n3]=size(X);

%% our method 
opts = [];
opts.maxIter=30;
opts.tol = 3e-3; % run to maxit by using negative tolerance %  1e-5     1e-3  3e-3
opts.Mtr = X; % pass the true tensor to calculate the fitting
opts.rank_min=round(5*ones(1,n3));opts.rank_min(1)=20;
EstCoreNway = round(20*ones(1,n3));EstCoreNway(1)=50;
Nway=[n1,n2,n3];
[~,~,TC,TotalRes,Time] =TCTFM_solver(data,known,Nway,EstCoreNway,opts,kk);
end

