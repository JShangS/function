function [ expA ] = fun_Expm( A )
%FUN_EXPM 此处显示有关此函数的摘要
%   此处显示详细说明
[UA, LA] = svd(A);
expA = UA * diag(exp(diag(LA))) * UA';
end

