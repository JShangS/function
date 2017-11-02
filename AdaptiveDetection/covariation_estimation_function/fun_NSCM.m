function [ R_NS ] = fun_NSCM( X )
%X:训练样本
%%归一化采样协方差矩阵
%一列是一个距离单元
[M,N] = size(X);
NX = zeros(M,N);
for i = 1:N
    NX(:,i) = X(:,i)/sqrt(X(:,i)'*X(:,i)/M);
end
R_NS = NX * NX'/L;
end

