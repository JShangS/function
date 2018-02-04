function [ Tanmf ] = fun_ANMF( R,x0,p )%%ACE
%%%Adaptive Nominalised Matched Filter
%%%R£ºCovariance Matrix£¬x0£ºCUT£¬p£ºsteering vector
%%%%%<<1995,Asymptotically Optimum Radar Detection in Compound-Gaussian Clutter>>
Tamf =  fun_AMF(R,x0,p);
iR = inv(R);
tmp=abs(x0'*iR*x0);
Tanmf = Tamf/tmp;
end

