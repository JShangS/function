function [ X ] = fun_TrainData_gauss(N,L,M)
%FUN_TRAINDATA 此处显示有关此函数的摘要
%   此处显示详细说明
%JerryShang，2017.11.09
%%产生训练数据
%%N：导向矢量维数
%%L：训练数据长度
if nargin<3
    error('输入参数至少为3个')
end
%%产生高斯杂波
X = (randn(N,L)+1i*randn(N,L))/sqrt(2);
M_half = M^0.5;
X = M_half*X;
end

