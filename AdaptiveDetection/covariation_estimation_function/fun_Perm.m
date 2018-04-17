function [ R ] = fun_Perm( X )
%FUN_PERM 此处显示有关此函数的摘要
%   此处显示详细说明
%%%%反对称结构协防差估计<Model Order Selection Rules for Covariance %
%%Structure Classification in Radar> 式（35）
[N,L]=size(X);
J = fun_permutation(N);
R = (X*X' + J*conj(X*X')*J)/2/L;
end

