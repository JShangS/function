function [ distance ] = fun_CholeskyDistance( A,B )
%FUN_REIMANDISTANCE 此处显示有关此函数的摘要
%   CholeskyDistance
%%% 根据文献<Covariance matrix estimation via geometric barycenters 
% and its application to radar training data selection>
distance = sqrt(trace((chol(A) - chol(B)) * (chol(A) - chol(B))'));
end

