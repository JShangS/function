function [ R_NS ] = fun_NSCM( X )
%%Adaptive matched filter detection in spherically invariant noise
%X:ѵ������
%%��һ������Э�������
%һ����һ�����뵥Ԫ
[M,N] = size(X);
NX = zeros(M,N);
for i = 1:N
    NX(:,i) = abs(X(:,i))/sqrt(abs(X(:,i)'*X(:,i)/M));
end
R_NS = abs(NX * NX'/N);
end

