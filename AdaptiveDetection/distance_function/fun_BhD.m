function [ distance ] = fun_BhD( A,B )
%FUN_BHD �˴���ʾ�йش˺�����ժҪ
%%%Geometric means and medians with applications to target detection(13)ʽ
%   Bhattacharyya distance
t1 = log(det((A + B)/2));
t2 = sqrt(det(A) * det(B));
distance = 2*sqrt(log(t1 / t2));
end

