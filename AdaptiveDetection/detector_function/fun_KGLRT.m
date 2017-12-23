function [ Tglrt ] = fun_KGLRT( R,x0,p )
%FUN_KGLRT 此处显示有关此函数的摘要
%   此处显示详细说明
%%Kelly的1S-GLRT
%%%R：协防差，x0：CUT，p：导向矢量
Tamf =  fun_AMF(R,x0,p);
iR = inv(R);
tmp=abs(x0'*iR*x0);
Tglrt = Tamf/(1+tmp);
end

