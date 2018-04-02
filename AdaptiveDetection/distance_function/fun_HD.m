function [ distance ] = fun_HD( A,B )
%FUN_BHD 此处显示有关此函数的摘要
%%%Geometric means and medians with applications to target detection(15)式
%   Hellinger distance
t1 = det(A)^(1/4) * det(B)^(1/4);
t2 = sqrt(det((A+B)/2));
distance = sqrt(2 - 2 * (t1 / t2));
end
