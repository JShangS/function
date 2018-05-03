function [ C ] = fun_corrcoef( X )
%FUN_CORRCOEF �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%���ϵ�����󣬼��㷽��Ϊ��Сǿ�����ģ�103ҳ
L = length(X);
c = zeros(1,L);
for i = 1:L % ��
    c(i) = sum(X(1:L-i+1) .* conj(X(1:L-i+1))) / (L-i+1);
end
C = toeplitz(c);
for i = 1:L
    for j = i+1:L
        C(i,j) = conj(C(i,j));
    end
end
end

