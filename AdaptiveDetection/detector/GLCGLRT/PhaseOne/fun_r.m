function [ r ] = fun_r( tau )
%FUN_R 此处显示有关此函数的摘要
%   此处显示详细说明
%%<Knowledge-Aided Bayesian Radar Detectors & Their Application to Live Data>
%%得到的自相关函数,phaseOne Radar
%%tau:自相关的差(i-j)
f0 = 1.23e9;
C=3e8;
lambda = C/f0;
w = 27/3.6;
d = 489.8*w^(-1.55)*f0^(-1.21);
ibeta = 0.1048*(log10(w)+0.4147);
r = d/(1+d) + 1/(1+d) *(1./(1+(4*pi*ibeta.*tau/lambda).^2));
end

