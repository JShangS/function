function [ T_MAM_GLRT ] = fun_MAM_GLRT( MAM,x0,p,lambda,mu )
%  文献 <Knowledge-based adaptive detection of radar targets in generalized Pareto clutter>
%%%的检测器，基于MAM模型，然后GLRT检测器在逆Gamma纹理下的最大似然估计
%%%%MAM:多模型
%%%p:目标导向矢量
%%%x0:CUT
%%%lambda：逆gamma纹理参数（形状参数）
%%%mu：逆gamma纹理参数（尺度参数）
if nargin<4 %%%%默认参数为3,1
    lambda = 3;
    mu = 1;
end
[N,~,L]=size(MAM); 
t1 = zeros(L,1);%%%1假设下系数
t0 = zeros(L,1);%%%0假设下系数
for i=1:L
    MAM(:,:,i) = inv(MAM(:,:,i));
    a=(p'*MAM(:,:,i)*x0)/(p'*MAM(:,:,i)*p);
    t1(i) = N/(lambda*mu*(x0-a*p)'*MAM(:,:,i)*(x0-a*p));
    t0(i) = N/(lambda*mu*x0'*MAM(:,:,i)*x0);
end
%%协方差估计
%%1假设下的
t1 = reshape(t1,1,1,L);
t1 = repmat(t1,N,N,1);
R_MAM1 = sum(t1.*MAM,3)/L;
%%0假设下的
t0 = reshape(t0,1,1,L);
t0 = repmat(t0,N,N,1);
R_MAM0 = sum(t0.*MAM,3)/L;
%%%检测器
a = (p'*R_MAM1*x0)/(p'*R_MAM1*p);
f1 = det(R_MAM1)*((x0-a*p)'*R_MAM1*(x0-a*p)+1/mu)^(-lambda-N);
f0 = det(R_MAM0)*(x0'*R_MAM0*x0+1/mu)^(-lambda-N);
T_MAM_GLRT = abs(f1/f0);
end

