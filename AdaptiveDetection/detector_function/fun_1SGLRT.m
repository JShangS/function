function [ T1sglrt ] = fun_1SGLRT( R,x0,p,mu )
%FUN_1SGRT 此处显示有关此函数的摘要
%   此处显示详细说明
%%<Radar detection based on compound-gaussian model with inverse gamma texture>
%%R:可是SCM，NSCM，CC等协方差
%%x0:CUT
%%p:steering vector
%%mu:invert gamma distribution,scale parameter.
Tamf =  fun_AMF(R,x0,p);
iR = inv(R);
tmp=abs(x0'*iR*x0);
T1sglrt = Tamf/(1/mu+tmp);
end

