function [ distance ] = fun_HD( A,B )
%FUN_BHD �˴���ʾ�йش˺�����ժҪ
%%%Geometric means and medians with applications to target detection(15)ʽ
%   Hellinger distance
t1 = det(A)^(1/4) * det(B)^(1/4);
t2 = sqrt(det((A+B)/2));
distance = sqrt(2 - 2 * (t1 / t2));
end
