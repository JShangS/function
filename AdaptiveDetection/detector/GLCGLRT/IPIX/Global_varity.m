clc
clear 
close all
%%%%%%%��������ȫ�ֱ�������ʹ��
cdfFile =  '19980205_170935_ANTSTEP.CDF';
cdfFile_t = cdfFile;
cdfFile_t(17:27)=[];
matFile = [cdfFile_t,'IPIX.mat'];