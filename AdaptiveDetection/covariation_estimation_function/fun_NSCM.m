function [ R_NS ] = fun_NSCM( X )
%X:ѵ������
%%��һ������Э�������
%һ����һ�����뵥Ԫ
[M,N] = size(X);
NX = zeros(M,N);
for i = 1:N
    NX(:,i) = X(:,i)/sqrt(X(:,i)'*X(:,i)/M);
end
R_NS = NX * NX'/L;
end

