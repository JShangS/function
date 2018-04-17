function [ value ] = fun_LogNorm( A )
%FUN_LOGNORM �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%%����LogDE�����ķ���
[M,N] = size(A);
%%%%%%%%%����A������������չ�ɾ���%%%%%%%%%
if M == 1 || N ==1
    A = A(:);
    A = A * A';
   [UA, LA] = svd(A);
   Evalue = abs(diag(LA));
   index_1 = find(Evalue<=1);
   Evalue(index_1) = 1;
   logmA = UA * diag(log(sqrt(Evalue))) * UA';
   value = abs(norm(logmA,'fro'));
else
   value = fun_LogED(A,eye(max(M,N)));
%     value = abs(norm(logm(A),'fro'));
end
% value = norm(logm(A),'fro');
end

