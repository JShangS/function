function [ R_NS ] = fun_NSCM( X )
%%Adaptive matched filter detection in spherically invariant noise
%X:训练样本
%%归一化采样协方差矩阵
%一列是一个距离单元
[M,N] = size(X);
NX = zeros(M,N);
for i = 1:N
    NX(:,i) = abs(X(:,i))/sqrt(abs(X(:,i)'*X(:,i)/M));
end
R_NS = abs(NX * NX'/N);
end

