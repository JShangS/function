clc 
clear
close all
load H067038_3iqHH_H067037_2iqVV.mat
N = 8;
PRF = 100;
Tr = 1/PRF;
Range = 30;
pulse = 1*N:2*N-1;
pulse_num = length(pulse);
Zhh_test = Zhh(pulse,:);
figure()
plot(abs(Zhh(30,:)))
figure(1)
imagesc(abs(Zhh))
[X,Y]=meshgrid(1:76,linspace(-1,1,pulse_num));
MTD = abs(fftshift(fft(Zhh_test,[],1)));
figure(2)
mesh(X,Y,MTD);
R_KA = zeros(N,N);
L_KA = 3800;
%%%%%%%%%%%%用CUT单元的所有数据计算先验协方差
for i = 1:L_KA
    X_KA = Zhh((i-1)*N+1:i*N,Range);
    R_KA = R_KA+X_KA*X_KA'/L_KA;
end
% X_KA = Zhh(40:end,Rang);
% L_KA = length(X_KA);
figure()
mesh(abs(R_KA))
% figure()
% R_test = Zhh_test(:,Range)*Zhh_test(:,Range)'/76;
% mesh(abs(R_test))

%%%用文献<Knowledge-Aided Bayesian Radar Detectors & Their Application to Live Data>
%%计算先验协方差
var_c = var(abs(Zhh(:,Range)));
for i=1:N
    for j=1:N
        r(i,j) = fun_r((i-j)*Tr);
    end
end
r = var_c*r;
R_KA=r;
figure()
mesh(abs(r))
%%%%%%%%%%%%%%%%%%%%%%



% R_KA = X_KA*X_KA'/L_KA;
save('Real_data_PhaseOneRadar_range30','R_KA','Zhh','Zhh_test','N','Range');
