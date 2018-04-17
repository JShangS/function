function [ R_LogCC,alpha ] = fun_LogCC( X,R_KA )
%FUN_LOGCC �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%X:ѵ������
%%R_KA:֪ʶЭ����
[M,N]=size(X);
R_SCMN = fun_SCMN(X);
rou = 0;
for i = 1:N
    Ri = X(:,i)*X(:,i)';
    Ri = fun_Positive(Ri);
    rou = rou + (fun_LogED( Ri , R_SCMN))^2/N^2;
end
alpha = rou/(rou + fun_LogED( R_SCMN , R_KA)^2);
R_LogCC = alpha * R_KA + (1 - alpha) * R_SCMN;

end
