function [ C ] = fun_corrcoef( X )
%FUN_CORRCOEF 此处显示有关此函数的摘要
%   此处显示详细说明
%%相关系数矩阵，计算方法为华小强大论文，103页
L = length(X);
c = zeros(1,L);
for i = 1:L % 行
    c(i) = sum(X(1:L-i+1) .* conj(X(1:L-i+1))) / (L-i+1);
end
C = toeplitz(c);
for i = 1:L
    for j = i+1:L
        C(i,j) = conj(C(i,j));
    end
end
end

