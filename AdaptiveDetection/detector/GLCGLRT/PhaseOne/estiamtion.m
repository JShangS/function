clc
clear 
close all
%%%%%%���Ʋ���
load H067038_3iqHH_H067037_2iqVV.mat
range = 72;
r = abs(Zvv(:,range));
[alpha,beta] = fun_IG_ML(r)
mesh(abs(Zvv))