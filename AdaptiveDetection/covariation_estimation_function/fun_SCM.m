function [ R_S ] = fun_SCM( X )
%采样协方差矩阵
%一列是一个距离单元
%X:训练样本
[~,N] = size(X);
R_S = (X*X');
end

