function [ R ] = fun_Perm( X )
%FUN_PERM �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%%%���ԳƽṹЭ�������<Model Order Selection Rules for Covariance %
%%Structure Classification in Radar> ʽ��35��
[N,L]=size(X);
J = fun_permutation(N);
R = (X*X' + J*conj(X*X')*J)/2/L;
end

