function [ Tclglrt ] = fun_iGLCGLRT( lambda,mu,RKA,R,x0,p )
%FUN_CLGRT 此处显示有关此函数的摘要
%   此处显示详细说明
%%色加载的GLRT的结果
[iR1,~] = fun_CLGLRT_icovariance(RKA,R,x0,p,1);%%H1下的协方差估计结果
[iR0,~] = fun_CLGLRT_icovariance(RKA,R,x0,p,0);%%H0下的协方差估计结果
[N,~]=size(x0);
a = (p'*iR1*x0)/(p'*iR1*p);
tmp1 = det(iR1)*((x0-a*p)'*iR1*(x0-a*p)+1/mu)^(-(N));%lambda
tmp0 = det(iR0)*(x0'*iR0*x0+1/mu)^(-(N));%lambda
Tclglrt = abs(tmp1/tmp0);
end

