function [ R_AML ] = fun_AML( X )
%%��Covariance matrixestimation for CFAR detection in correlated heavy tailed clutter��
%FUN_AML �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%���������ȻЭ����
%X:ѵ������
[M,N] = size(X);
R_AML = eye(M,M);%R_x0;%eye(N,N);%�Ե�λ��Ϊ������ֵ�ǣ��ڶ��ε������ΪNSCM���
tao_child = 1;%%���ε���ֵ
tao_parent = 0;%%�ϴε���ֵ
% count = 1;
while (abs(tao_child-tao_parent)>0.01)%
    tao_parent = tao_child;
    iR_AML = inv(R_AML);
    tao_child = diag(abs(X'*iR_AML*X)/M);
    R_AML_t = 0;
    for i = 1:N
        R_AML_t = R_AML_t+X(:,i)*X(:,i)'/N/tao_child(i);
    end
    R_AML = (R_AML_t);
end
end

