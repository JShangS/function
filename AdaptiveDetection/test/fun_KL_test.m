function [ R,alpha ] = fun_KL_test( R_KA,R1,R_estimation )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[N,~] = size(R_KA);
t1 = R_estimation * (R_KA - R1) - R1;
t2 = R_KA - R1;
alpha = trace(t1/t2)/N;
alpha = max(min(1,alpha),0);
R = alpha * R_KA + (1-alpha) * R1;
end

