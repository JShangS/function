function [ distance ] = fun_LogED( A,B )
%FUN_REIMANDISTANCE 此处显示有关此函数的摘要
%Log Euclidean Distance
%   此处显示详细说明
%%% 根据文献<Covariance matrix estimation via geometric barycenters 
% and its application to radar training data selection>
[UA, LA] = svd(A);
logA = UA * diag(log(diag(LA))) * UA';
[UB, LB] = svd(B);%%特征向量，特征值，
logB = UB * diag(log(diag(LB))) * UB';
distance = (trace((logA - logB) * (logA - logB)'));
% distance = (trace((logm(A) - logm(B)) * (logm(A) - logm(B))'));
end

