function [ J ] = fun_permutation( N )
%FUN_PERMUTATION �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%%�����û����󣬷��Խ��߶�Ϊ1����Ϊ0
J = zeros(N,N);
for i = 1:N
    J(i,N-i+1) = 1;
end
end

