function [ a ] = fun_test( R0,R1,Rcut )
%FUN_TEST 此处显示有关此函数的摘要
%   此处显示详细说明
[N,~] = size(R0);
a = 0;
Rcut = fun_Positive(Rcut);
iter_num = 100;
for i = 1:iter_num
    t1 = det((a*R0 + (1-a)*R1) \ Rcut);
    a =  abs(trace((Rcut*t1 - R1)/(R0 - R1))/N)
end
% a = min(1,a);
end

