function [ Th_AMF ] = fun_Th_AMF( K,N,PFA )
%FUN_TH_ACE 此处显示有关此函数的摘要
%   此处显示详细说明
%K：训练样本列数
%N：检测单元行数
%PFA:虚警率
Th_AMF = (K+1)/(K-N+1)*((PFA)^(1/(K-N+2))-1)
end

