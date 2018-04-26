function [ R_AML, alpha ] = fun_AMLCC2( X,R_KA )
%%��Covariance matrixestimation for CFAR detection in correlated heavy tailed clutter��
%FUN_AML �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%���������ȻЭ����
%X:ѵ������
[M,N] = size(X);
R_AML = eye(M,M);%R_x0;%eye(N,N);%�Ե�λ��Ϊ������ֵ�ǣ��ڶ��ε������ΪNSCM���
R_AML0 = zeros(M,M);
% count = 1;
alpha0 = 1;
alpha = 0;
while (abs(alpha0-alpha)>1e-2)%
%     R_AML_inv = inv(R_AML);
    alpha0=alpha;
    R_AML0 = R_AML;
    tao_child = diag((X'/R_AML0*X)/M);
    t1 = 0;
    t2 = 0;
    t3 = 0;
    for i = 1:N
        t1 = t1 + norm(X(:,i) * X(:,i)'/tao_child(i)-R_AML0,'fro')^2;
        t2 = t2 + norm(R_KA - X(:,i) * X(:,i)'/tao_child(i),'fro')^2;
        t3 = t3 + X(:,i) * X(:,i)'/N/tao_child(i);
    end
    alpha = max(min(1,t1/t2),0);
    R_AML = alpha * R_KA + (1-alpha) * t3;
%     R_AML = t3;
end
end

