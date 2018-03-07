function [ distance ] = fun_RootED( A, B )
%FUN_REIMANDISTANCE 此处显示有关此函数的摘要
%Root-Euclidean distance
%   此处显示详细说明
%%% 根据文献<Covariance matrix estimation via geometric barycenters 
% and its application to radar training data selection>
[UA, LA] = eig(A);
A2 = UA * LA^0.5 * UA';
[UB, LB] = eig(B);
B2 = UB * LB^0.5 * UB';
distance = sqrt(trace((A2-B2) * (A2-B2)'));
end

