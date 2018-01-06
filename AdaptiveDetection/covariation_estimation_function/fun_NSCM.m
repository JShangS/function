function [ R_NS ] = fun_NSCM( X,opt )
%X:训练样本
%%归一化采样协方差矩阵
%一列是一个距离单元
if nargin==1
    opt=1;
end
[M,N] = size(X);
NX = zeros(M,N);
R_NS = zeros(M,M);
if opt==1
    %%Adaptive matched filter detection in spherically invariant noise
    %%数据归一化
    for i = 1:N
        NX(:,i) = X(:,i)/sqrt(norm(X(:,i),'fro')^2/M);
    end
    R_NS = (NX * NX');
elseif opt==2
    %%Performance analysis of two covariance matrix estimators in
    %%compound-Gaussian clutter
    for i = 1:N
         t = X(:,i)*X(:,i)'/(X(:,i)'*X(:,i))/M;
        R_NS = R_NS + t;
    end   
end
end

