function [ PSD ] = fun_PSD( R,ft,N )
% Capon PSD estimator
%%《 Estimation of the Covariance Matrix Based on Multiple A-Priori Models, A. De Maio》
%%N：归一化多普勒点数
%%R：协方差矩阵
[M,~] = size(R);
if nargin == 1
    N = 1000;
    ft = linspace(-0.5,0.5,N);
elseif nargin == 2
    N = length(ft);
else
    error('参数至少一个协方差')
end
nn = (0:M-1)';
PSD = zeros(N,1);
iR = inv(R);
for i = 1:N
    p = exp(1j*2*pi*nn*ft(i))/sqrt(M);
    PSD(i) = (1/(p'*iR*p));
end
PSD=abs(PSD)/max(abs(PSD));
% PSD = 10*log10(PSD);
PSD=(fun_value2dB(PSD));
end

