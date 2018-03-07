function [ distance ] = fun_RootED( A, B )
%FUN_REIMANDISTANCE �˴���ʾ�йش˺�����ժҪ
%Root-Euclidean distance
%   �˴���ʾ��ϸ˵��
%%% ��������<Covariance matrix estimation via geometric barycenters 
% and its application to radar training data selection>
[UA, LA] = eig(A);
A2 = UA * LA^0.5 * UA';
[UB, LB] = eig(B);
B2 = UB * LB^0.5 * UB';
distance = sqrt(trace((A2-B2) * (A2-B2)'));
end

