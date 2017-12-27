function [ RW ] = fun_Wishart_R(R,mu)
%FUN_WISHART_R 此处显示有关此函数的摘要
%   此处显示详细说明
%%%服从均值为R的wishart分布的协方差。，R可由fun_rho得到，也可以是先验的协方差。
%%<Conjugate Bayesian analysis of the Gaussian distribution>
%%%freedom自由度
RW = wishrnd(1/mu*R,mu);
end

