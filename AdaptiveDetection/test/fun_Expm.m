function [ expA ] = fun_Expm( A )
%FUN_EXPM �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[UA, LA] = svd(A);
expA = UA * diag(exp(diag(LA))) * UA';
end

