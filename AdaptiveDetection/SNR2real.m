function [ y ] = SNR2real( SNR,scaler )
%SNR2REAL 此处显示有关此函数的摘要
%   SNR换算成正常值
y = 10.^(SNR./scaler);
end

