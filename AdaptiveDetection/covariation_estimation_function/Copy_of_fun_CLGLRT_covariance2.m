function [ R_CLGRT,alpha,beta ] = fun_CLGLRT_covariance2(RKA,R,x0,p,opt )
%FUN_CLGLRT_CON 此处显示有关此函数的摘要
%   此处显示详细说明
%lamda:逆gamma分布的形状参数，mu：尺度参数
%RKA：先验协方差，R:采样协方差
%x0：待检测单元
%opt：H1，H0假设下的协方差估计结果
%p：导向矢量
iRKA = inv(RKA);
iR = inv(R);
[N,~] = size(x0);
if opt == 1 %%H1假设下的协方差估计
   a1 = (p'*iRKA*x0)/(p'*iRKA*p);
   a2 = (p'*iR*x0)/(p'*iR*p);
   fai1 = (((x0-a1*p)'*iRKA*(x0-a1*p))/((x0-a2*p)'*iR*(x0-a2*p)));%fai1 = (((x0-a1*p)'*iRKA*(x0-a1*p))/((x0-a2*p)'*iR*(x0-a2*p)));
   alpha1 = (1/(1+fai1));
   beta1 = (fai1/(1+fai1));
   %%%R_CLGRT
   R_CLGRT = alpha1*RKA+beta1*R;
   alpha = alpha1;
   beta = beta1;
%     R_CLGRT = alpha1*RKA + beta1*R;
elseif opt == 0%H0假设下的协方差估计
    fai0 =  ((x0'*iRKA*x0)/(x0'*iR*x0));%fai0 =  ((x0'*iRKA*x0)/(x0'*iR*x0));
    alpha0 = (1/(1+fai0));
    beta0 = (fai0/(1+fai0));
    R_CLGRT = alpha0*RKA+beta0*R;
    alpha = alpha0;
    beta = beta0;
end
end

