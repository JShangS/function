function [ distance ] = fun_EuclideanDistance( A,B )
%FUN_REIMANDISTANCE 此处显示有关此函数的摘要
%   此处显示详细说明
%%% 根据文献<Covariance matrix estimation via geometric barycenters 
% and its application to radar training data selection>
distance = sqrt(trace((A-B) * (A-B)'));
end

