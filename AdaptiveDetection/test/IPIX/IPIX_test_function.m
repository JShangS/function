clc
clear 
close all
n = 2;
SNRout=-5:1:25; % Êä³öSNR
PFA=1e-2;% PFA=1e-4;

f{1,1} = @(X)fun_SCMN(X);
f{1,2} = @(R,x0,s)fun_ANMF(R,x0,s);
f{1,3} = 'ANMF with SCM';

% f{2,1} = @(X)fun_SCMN(X);
% f{2,2} = @(R,x0,s)fun_AMF(R,x0,s);
% f{2,3} = 'AMF with SCM';

[Pd, Th] = fun_IPIX_Detction(f,n,SNRout,PFA);