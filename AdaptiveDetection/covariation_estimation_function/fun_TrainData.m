function [ X ] = fun_TrainData(N,L, M,opt)
%FUN_TRAINDATA 此处显示有关此函数的摘要
%   此处显示详细说明
%%产生训练数据
%%N：导向矢量维数
%%L：训练数据长度
%%opt：产生的训练数据类型：1：高斯；2：未完待续....
if nargin==3
    opt = 1;
end
if nargin<3
    error('输入参数至少为3个')
end
%%产生高斯杂波
if opt == 1
    X = (randn(N,L)+1i*randn(N,L))/sqrt(2);
    M_half = M^0.5;
    X = M_half*X;
end
end

