function [ Train ] = fun_TrainData_K( N,L, M, v)
%FUN_TRAINDATA_K 此处显示有关此函数的摘要
%   此处显示详细说明
%%K纹理分布的复合高斯分布
%%N,L行列
%%M，复高斯协方差
[U,D] = eig(M);
M_half = U*D.^0.5;
Train=zeros(N,L);
for l=1:L
    gs=sqrt(1/2)*(randn(N,1)+1j*randn(N,1)); 
    taos=gamrnd(v,1/v,1,1);  %K-distribution clutter
    Train(:,l)=sqrt(taos)*(M_half*gs);
end
end

