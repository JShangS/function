function [ distance ] = fun_LogED( A,B )
%FUN_REIMANDISTANCE �˴���ʾ�йش˺�����ժҪ
%Log Euclidean Distance
%   �˴���ʾ��ϸ˵��
%%% ��������<Covariance matrix estimation via geometric barycenters 
% and its application to radar training data selection>
[UA, LA] = svd(A);
logA = UA * diag(log(diag(LA))) * UA';
[UB, LB] = svd(B);%%��������������ֵ��
logB = UB * diag(log(diag(LB))) * UB';
distance = (trace((logA - logB) * (logA - logB)'));
% distance = (trace((logm(A) - logm(B)) * (logm(A) - logm(B))'));
end

