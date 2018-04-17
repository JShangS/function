function [ value ] = fun_LogNorm( A )
%FUN_LOGNORM 此处显示有关此函数的摘要
%   此处显示详细说明
%%%基于LogDE导出的范数
[M,N] = size(A);
%%%%%%%%%若果A是列向量则扩展成矩阵%%%%%%%%%
if M == 1 || N ==1
    A = A(:);
    A = A * A';
   [UA, LA] = svd(A);
   Evalue = abs(diag(LA));
   index_1 = find(Evalue<=1);
   Evalue(index_1) = 1;
   logmA = UA * diag(log(sqrt(Evalue))) * UA';
   value = abs(norm(logmA,'fro'));
else
   value = fun_LogED(A,eye(max(M,N)));
%     value = abs(norm(logm(A),'fro'));
end
% value = norm(logm(A),'fro');
end

