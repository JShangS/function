function [ y ] = SNR2real( SNR,scaler )
%SNR2REAL �˴���ʾ�йش˺�����ժҪ
%   SNR���������ֵ
y = 10.^(SNR./scaler);
end

