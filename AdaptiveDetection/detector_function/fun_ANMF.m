function [ Tanmf ] = fun_ANMF( R,x0,p )%%ACE
%%%Adaptive Nominalised Matched Filter
%%%R��Covariance Matrix��x0��CUT��p��steering vector
%%%%%<<1995,Asymptotically Optimum Radar Detection in Compound-Gaussian Clutter>>
Tamf =  fun_AMF(R,x0,p);
tmp=abs(x0'*iR*x0);
Tanmf = Tamf/tmp;
end

