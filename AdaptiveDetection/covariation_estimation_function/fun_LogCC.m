function [ R_CC,alpha0 ] = fun_LogCC(X,R,R_KA )
%FUN_CC �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%������⣬ѵ���������Ƶ�Э���������Э�����������ϣ�����͹�Ż��õ����ϵ����
%%X:ѵ������
%R,�������Ƶ�Э����
%R_KA:����Э����
% R = fun_SCMN(X);
R_CC = R;
[N,K] = size(X);
Rim = zeros(N,N,K);

for i = 1:K
    Rim(:,:,i) = logm(fun_Positive( X(:,i),4));
    
end
alpha = 0;
alpha0 =0;
for k = 1:20
    alpha = alpha0;
    t1 = 0;
    t2 = 0;
%     t3 = zeros(N,N);
    for i = 1:K 
        t1 = t1 + norm( Rim(:,:,i) - logm(R_CC),'fro')^2;
        t2 = t2 + norm(logm(R_KA) - Rim(:,:,i),'fro')^2;
    end
    alpha0 = max(min(1,(t1/t2)),0);
%     for i = 1:K
%         t3 = t3 + alpha0 * R_KA +(1-alpha0) * X(:,i) * X(:,i)';
%     end
%     R_CC_0 = R_CC;
%     R_CC = t3/K;    
    R_CC_0 = R_CC;
    R_CC = expm(alpha0 * logm(R_KA) +(1-alpha0) * logm(R_CC_0));
    if abs(alpha0- alpha)<1e-2
        break;
    end
end
end


