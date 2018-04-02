function [ distance ] = fun_BhD( A,B )
%FUN_BHD 此处显示有关此函数的摘要
%%%Geometric means and medians with applications to target detection(13)式
%   Bhattacharyya distance
t1 = log(det((A + B)/2));
t2 = sqrt(det(A) * det(B));
distance = 2*sqrt(log(t1 / t2));
end

