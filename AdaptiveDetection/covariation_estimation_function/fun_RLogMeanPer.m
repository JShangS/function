function [ R ] = fun_RLogMeanPer( X )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%%Log-E�ľ�ֵЭ����,���ڷ��Գ�����
[N,L] = size(X);
J = fun_permutation(N);
logmR = 0;
for i = 1:L
%     Ri = X(:,i) * X(:,i)';
    Ri = fun_Positive(X(:,i),3);
    logmR = logmR + fun_Logm(Ri + J * conj(Ri) * J);
end
% logmR = logm(X*X' + J * conj(X*X') * J) + logm(0.5*eye(N));

R = fun_Expm(logmR/L+logm(0.5*eye(N)));
end


