function DT= DTRTC_Video( X,Omega,Nway2 )

%% produce data
data=X(Omega);known=Omega;
Nway{1}=size(X);TC{1}=zeros(Nway{1});TC{1}(known)=data;
TC{2}=reshape(Unfold(TC{1},Nway{1},3),Nway2);
Nway{2}=size(TC{2});
%% our method 
opts = [];
opts.maxIter=150;
opts.tol = 1e-4; % run to maxit by using negative tolerance
opts.Mtr = X; % pass the true tensor to calculate the fitting
opts.TC = TC;
opts.gamma=10;
opts.rank_min{1} = 1*ones(1,Nway{1}(3));opts.rank_min{1}(1)=50; %50
opts.rank_min{2} = 1*ones(1,Nway{2}(3));opts.rank_min{2}(1)=10; %10
EstCoreNway{1} = 70*ones(1,Nway{1}(3));EstCoreNway{1}(1) =120;
EstCoreNway{2} = 10*ones(1,Nway{2}(3));   %10
[~,~,DT] = DTRTC_solver(data,known,Nway,EstCoreNway,opts);
end

