function [ J ] = fun_permutation( N )
%FUN_PERMUTATION 此处显示有关此函数的摘要
%   此处显示详细说明
%%%产生置换矩阵，反对角线都为1其余为0
J = zeros(N,N);
for i = 1:N
    J(i,N-i+1) = 1;
end
end

