function DT= DTRTC_Color( X,Omega,Nway2 )

%% produce data
data=X(Omega);known=Omega;
Nway{1}=size(X);TC{1}=randn(Nway{1});TC{1}(known)=data;
TC{2}=reshape(Unfold(TC{1},Nway{1},3),Nway2); 
Nway{2}=size(TC{2});
%% our method 
opts = [];
opts.maxIter=100;
opts.tol = 5e-4; % run to maxit by using negative tolerance
opts.Mtr = X; % pass the true tensor to calculate the fitting
opts.TC = TC;
opts.gamma=0.1;
opts.alpha_adj = 0;
opts.rank_adj = -1*ones(1,Nway{1}(3));
opts.rank_inc = 1*ones(1,Nway{1}(3));opts.rank_max = [50,20,20];
%% picture
opts.rank_min{1} = 5*ones(1,Nway{1}(3));opts.rank_min{1}(1)=30;
opts.rank_min{2} = 1*ones(1,Nway{2}(3));
EstCoreNway{1} = 30*ones(1,Nway{1}(3));EstCoreNway{1}(1)=200;
EstCoreNway{2} = 3*ones(1,Nway{2}(3));
[~,~,DT] = DTRTC_solver(data,known,Nway,EstCoreNway,opts);
end

