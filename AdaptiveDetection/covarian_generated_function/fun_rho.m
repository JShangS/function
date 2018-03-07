function [ R_rho ] = fun_rho( rho,N,fd)
%FUN_ROU 此处显示有关此函数的摘要
%   此处显示详细说明
%%最基本的协方产生。
%%%rho，迟滞因子
%%%N，导向矢量维数
R_rho = zeros(N,N);
if nargin<3
    fd = 0;
end
L = length(rho);
for l = 1:L
    for i=1:N
        for j=1:N
            R_rho(i,j,l)=rho(l)^abs(i-j)*exp(1j*2*pi*fd*(i-j));
        end
    end
end
end


