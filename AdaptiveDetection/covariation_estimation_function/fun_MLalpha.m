function [ R_MLalpha, Alpha] = fun_MLalpha( X,R_KA,x0)
%��Knowledge-Aided Space-Time Adaptive Processing,2011��
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%ɫ����ϵ��alpha�������Ȼ����,
%%X��ѵ������
%%R:ѵ������Э����
%%R_KA:����Э����
%%x0��Ŀ���ⵥԪ
[M,N] = size(X);
R = fun_SCMN(X);
G = X'*inv(R_KA)*X;
[U,A] = eig(G);%%U:�Ͼ���A�Խ���,G=U*A*U';
P = inv(R_KA)*X*U;
R_KA_inv = inv(R_KA);
alpha_max = 0;
C_max = 0; 
X=cat(2,X,x0);
for alpha = 0.01:0.01:1
    eta = (1-alpha)/alpha/N;
    C_t1 = log(alpha^M*det(R_KA)*det(eye(N,N)+eta*A));
    C_t2 = x0'*(1/alpha*(R_KA_inv-eta*P*inv(eye(N,N)+eta*A)*P'))*x0;
    C_t = -(C_t1+C_t2);
%     C_t = 0;
%     for i = 1:N+1
%         C_t2 =X(:,i)'*(1/alpha*(R_KA_inv-eta*P*inv(eye(N,N)+eta*A)*P'))*X(:,i);
%         C_t =  C_t+(C_t1+C_t2);
%     end
%     C_t = -C_t;


%     C_t1 = log(det(alpha*R_KA+(1-alpha)*R_KA));
%     C_t2 =x0'*(alpha*R_KA+(1-alpha)*R_KA)*x0;
%     C_t =  -(C_t1+C_t2);
    if (C_t)>(C_max)
        alpha_max = alpha;
        C_max = C_t;
    end
end
Alpha = alpha_max;
R_MLalpha = Alpha*R_KA + (1-Alpha)*R;
end

