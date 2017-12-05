function [ R_S ] = fun_SCMN( X )
%采样协方差矩阵/N
%一列是一个距离单元
%X:训练样本
[~,N] = size(X);
R_S = (X*X'/N);
end

