function [ distance ] = fun_EuclideanDistance( A,B )
%FUN_REIMANDISTANCE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%% ��������<Covariance matrix estimation via geometric barycenters 
% and its application to radar training data selection>
distance = sqrt(trace((A-B) * (A-B)'));
end

