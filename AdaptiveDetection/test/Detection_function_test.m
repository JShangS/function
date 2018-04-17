clc
clear 
close all
n = 2;
SNRout=-5:1:25; %  ‰≥ˆSNR
PFA=1e-1;% PFA=1e-4;
rou = 0.95;
N = 8;
sigma_t = 0.1;


rouR = fun_rho(0.95,N);
R_KA = zeros(size(rouR));
for i = 1:1000
    t = normrnd(1,sigma_t,N,1);%%0~0.5%% ß≈‰œÚ¡ø
    R_KA = R_KA + rouR.*(t*t')/1000;
end

f{1,1} = @(X)fun_SCMN(X);
f{1,2} = @(R,x0,s)fun_ANMF(R,x0,s);
f{1,3} = 'ANMF with SCM';

f{2,1} = @(X)fun_CC(X,R_KA);
f{2,2} = @(R,x0,s)fun_ANMF(R,x0,s);
f{2,3} = 'ANMF with CC';

f{3,1} = @(X)fun_CCIter(X,R_KA);
f{3,2} = @(R,x0,s)fun_ANMF(R,x0,s);
f{3,3} = 'ANMF with CCIter';

[M,L] = size(f);
[Pd, Th] = fun_Detection(f,n,SNRout,PFA,N,rouR);

figure(1)
hold on
for i = 1:M
    plot(Pd(i,:))
end
