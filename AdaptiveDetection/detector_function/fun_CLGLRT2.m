function [ Tclglrt ] = fun_CLGLRT2( lambda,mu,RKA,R,x0,p )
%FUN_CLGRT 此处显示有关此函数的摘要
%   此处显示详细说明
%%色加载的GLRT的结果
R1 = fun_CLGLRT_covariance(lambda,mu,RKA,R,x0,p,1);%%H1下的协方差估计结果
R0 = fun_CLGLRT_covariance(lambda,mu,RKA,R,x0,p,0);%%H0下的协方差估计结果
[N,~]=size(x0);
iR1 = inv(R1);
iR0 = inv(R0);
a = (p'*iR1*x0)/(p'*iR1*p);
tmp1 = det(R0)*(x0'*iR0*x0)^(N); %,;det(R0)*((x0-a*p)'*iR1*(x0-a*p)+1/mu)^(-(lambda+N));
tmp0 = det(R1)*((x0-a*p)'*iR1*(x0-a*p))^(N);%,;det(R1)*(x0'*iR0*x0+1/mu)^(-(lambda+N)); 
Tclglrt = abs(tmp1/tmp0);
end

