function [ y ] = Qinv( x )
%QINV 此处显示有关此函数的摘要
%   此处显示详细说明
%  右尾函数反函数，x为概率，右尾标准正态分布的积分下限P(Y>y)=x
y = sqrt(2) * erfinv(1-2*x);
end

