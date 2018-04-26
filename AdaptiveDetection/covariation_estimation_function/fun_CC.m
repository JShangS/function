function [ R_CC,alpha0] = fun_CC( X,R,R_KA )
%FUN_CC 此处显示有关此函数的摘要
%《On Using a priori Knowledge in Space-Time Adaptive Processing》
%   此处显示详细说明
%%训练样本估计的协方差和先验协方差的线性组合，利用凸优化得到组合系数。
%%X:训练样本
%R,样本估计的协方差
% R = abs(fun_SCMN(X));
%R_KA:先验协方差
[M,N]=size(X);
rou_ba_t = 0;
for i = 1:N
    rou_ba_t = rou_ba_t+norm(X(:,i),'fro')^4/(N^2);
end
rou_ba = rou_ba_t-norm(R,'fro')^2/N;%R_KA
alpha0 = rou_ba/(rou_ba+norm(R-R_KA,'fro')^2);
alpha0 = max(min(1,alpha0),0);
% rou_ba = sum(diag(X'*X).^2)/N^2-sum(sum(abs(R).^2))/N;%（18）式,
% alpha0 = rou_ba/(rou_ba+sum(sum(abs(R-R_KA).^2)));
R_CC = (1-alpha0)*R+alpha0*R_KA;
end

