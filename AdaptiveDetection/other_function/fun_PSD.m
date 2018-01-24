function [ PSD ] = fun_PSD( R,ft,N )
%PSD �˴���ʾ�йش˺�����ժҪ
% Capon PSD estimator
%%�� Estimation of the Covariance Matrix Based on Multiple A-Priori Models, A. De Maio��
%%N����һ�������յ���
%%R��Э�������
[M,~] = size(R);
if nargin == 1
    N = 1000;
    ft = linspace(-0.5,0.5,N);
elseif nargin == 2
    N = length(ft);
else
    error('��������һ��Э����')
end
nn = (0:M-1)';
PSD = zeros(N,1);
iR = inv(R);
for i = 1:N
    p = exp(1j*2*pi*nn*ft(i));
    PSD(i) = 1/(p'*iR*p);
end
end

