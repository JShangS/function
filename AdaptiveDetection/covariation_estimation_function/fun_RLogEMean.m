function [ R ] = fun_RLogEMean( X )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%%Log-E�ľ�ֵЭ����
[N,L] = size(X);
logm_R = zeros(N,N);
for i = 1:L
%     Ri = X(:,i) * X(:,i)';
    Ri = fun_Positive(X(:,i),3);
    logm_R = logm_R + fun_Logm(Ri);
end
logm_R = logm_R/L;
R = fun_Expm(logm_R);
end

