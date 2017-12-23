function [ Tamf ] = fun_AMF( R,x0,p )
%FUN_AMF 此处显示有关此函数的摘要
%   此处显示详细说明
%%AMF
%%%R：协防差，x0：CUT，p：导向矢量
iR = inv(R);
Tamf = abs(p'*iR*x0)^2/abs(p'*iR*p);   
end

