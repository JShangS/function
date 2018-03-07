function [ distance ] = fun_LogED( A,B )
%FUN_REIMANDISTANCE 此处显示有关此函数的摘要
%Log Euclidean Distance
%   此处显示详细说明
%%% 根据文献<Covariance matrix estimation via geometric barycenters 
% and its application to radar training data selection>
[UA, LA] = eig(A);
logA = UA * log(LA) * UA';
[UB, LB] = eig(B);
logB = UB * log(LB) * UB';
distance = sqrt(trace((logA - logB) * (logA - logB)'));
end

