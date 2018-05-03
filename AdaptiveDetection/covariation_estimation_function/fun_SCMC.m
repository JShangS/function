function [ R_S ] = fun_SCMC( X  )
%相关系数矩阵/N
%一列是一个距离单元
%X:训练样本
[M,N] = size(X);
R_S = zeros(M,M);
for i = 1:N
    R_S = R_S + fun_Positive(X(:,i),4)/N;
end
end

