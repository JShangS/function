function [ Train ] = fun_TrainData_IGCC( N,L,M,lamda,mu )
%FUN_TRAINDATA_IGCC 此处显示有关此函数的摘要
%   此处显示详细说明
%%逆gamma纹理的复合高斯杂波
%%N,L行列,导向矢量维度，参考单元数目
%%M，复高斯协方差
%%lamda：形状参数
%%mu：尺度参数
pct=sqrt(ones(L,1)./gamrnd(lamda,mu,L,1));
X = fun_TrainData_gauss(N,L,M);
Train=(ones(N,1)*pct.').*X;
end

