function [ a ] = fun_test( R0,R1,Rcut )
%FUN_TEST �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
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

