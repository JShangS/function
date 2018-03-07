function [ Topt ] = fun_OPT( R,x0,p )
%%%最优检测器，协方差矩阵已知的情况
%%%R：Covariance Matrix，x0：CUT，p：steering vector
iR = inv(R);
Topt = abs(p'*iR*x0)^2 / (p'*iR*p) / (x0'*iR*x0);
end

