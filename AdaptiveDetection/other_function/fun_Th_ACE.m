function [ Th_ACE ] = fun_Th_ACE( K,N,PFA )
%FUN_TH_ACE 此处显示有关此函数的摘要
%   此处显示详细说明
%K：训练样本列数
%N：检测单元行数
%PFA:虚警率
t1 = 1-(PFA)^(1/(K-N+1));
t2 = 1-(K-N+1)/(K+1)*(PFA)^(1/(K-N+1));
Th_ACE = t1/t2;
end

