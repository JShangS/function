function [ R_logAML ] = fun_LogAML( X )
%FUN_LOGAML �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[M,N] = size(X);
R_logAML = eye(M,M);%R_x0;%eye(N,N);%�Ե�λ��Ϊ������ֵ�ǣ��ڶ��ε������ΪNSCM���
tao_child = 1;%%���ε���ֵ
tao_parent = 0;%%�ϴε���ֵ
% count = 1;
while (abs(tao_child-tao_parent)>0.5)%
    tao_parent = tao_child;
    iR_AML = inv(R_logAML);
    tao_child = diag((X'*iR_AML*X)/M);
    R_AML_t = 0;
    for i = 1:N
%         tao = (X(:,i)'*R_AML_inv*X(:,i))/M;
        Ri = fun_Positive(X(:,i),4);
        R_AML_t = R_AML_t+logm(Ri/tao_child(i));
    end
    R_AML_t = expm(R_AML_t/N);
%     if norm(R_AML_t-R_logAML,'fro')<1e-2
%         break;
%     end
    R_logAML = R_AML_t;
end

end

