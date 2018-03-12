function [ distance ] = fun_LogED( A,B )
%FUN_REIMANDISTANCE �˴���ʾ�йش˺�����ժҪ
%Log Euclidean Distance
%   �˴���ʾ��ϸ˵��
%%% ��������<Covariance matrix estimation via geometric barycenters 
% and its application to radar training data selection>
[UA, LA] = eig(A);
logA = UA * diag(log(diag(LA))) * UA';
[UB, LB] = eig(B);%%��������������ֵ��
logB = UB * diag(log(diag(LB))) * UB';
distance = (trace((logA - logB) * (logA - logB)'));
end

