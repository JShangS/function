function [ R_S ] = fun_SCMC( X  )
%���ϵ������/N
%һ����һ�����뵥Ԫ
%X:ѵ������
[M,N] = size(X);
R_S = zeros(M,M);
for i = 1:N
    R_S = R_S + fun_Positive(X(:,i),4)/N;
end
end

