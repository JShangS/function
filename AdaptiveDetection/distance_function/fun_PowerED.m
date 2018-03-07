function [ distance ] = fun_PowerED( A, B, alpha )
%FUN_REIMANDISTANCE �˴���ʾ�йش˺�����ժҪ
%Power-Euclidean distance
%   �˴���ʾ��ϸ˵��
%%% ��������<Covariance matrix estimation via geometric barycenters 
% and its application to radar training data selection>
if nargin<3
    alpha = 2;
end
[UA, LA] = eig(A);
A2 = UA * LA^alpha * UA';
[UB, LB] = eig(B);
B2 = UB * LB^alpha * UB';
distance = sqrt(trace((A2-B2) * (A2-B2)'));
end

