function [ R_NS ] = fun_NSCMN( X )
%%Adaptive matched filter detection in spherically invariant noise
%X:训练样本
%%归一化采样协方差矩阵/N
%一列是一个距离单元
[M,N] = size(X);
NX = zeros(M,N);
%%数据归一化
for i = 1:N
    NX(:,i) = X(:,i)/sqrt(norm(X(:,i),'fro')^2/M);
end
R_NS = (NX * NX'/N);
end

