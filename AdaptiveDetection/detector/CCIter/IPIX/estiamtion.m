clc
clear 
close all
%%%%%%¹À¼Æ²ÎÊý
Read_Display_Data
range = 14;
% sig = sig.';
r = sig(:,range);
[alpha,beta] = fun_IG_ML(r)
mesh(abs(sig))