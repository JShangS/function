function [ distance ] = fun_CholeskyDistance( A,B )
%FUN_REIMANDISTANCE �˴���ʾ�йش˺�����ժҪ
%   CholeskyDistance
%%% ��������<Covariance matrix estimation via geometric barycenters 
% and its application to radar training data selection>
distance = sqrt(trace((chol(A) - chol(B)) * (chol(A) - chol(B))'));
end

