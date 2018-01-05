clc
clear 
close all
%%%%%%%声明几个全局变量方便使用
cdfFile =  '19980223_170435_ANTSTEP.CDF';
cdfFile_t = cdfFile;
cdfFile_t(17:27)=[];
matFile = [cdfFile_t,'IPIX.mat'];