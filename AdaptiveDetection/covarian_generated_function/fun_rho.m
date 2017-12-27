function [ R_rho ] = fun_rho( rho,N)
%FUN_ROU 此处显示有关此函数的摘要
%   此处显示详细说明
%%最基本的协方产生。
%%%rho，迟滞因子
%%%N，导向矢量维数
R_rho = zeros(N,N);
for i=1:N
    for j=1:N
        R_rho(i,j)=rho^abs(i-j);
    end
end
end

