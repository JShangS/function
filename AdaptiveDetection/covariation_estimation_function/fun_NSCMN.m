function [ R_NS ] = fun_NSCMN( X )
%%Adaptive matched filter detection in spherically invariant noise
%X:ѵ������
%%��һ������Э�������/N
%һ����һ�����뵥Ԫ
[M,N] = size(X);
NX = zeros(M,N);
%%���ݹ�һ��
for i = 1:N
    NX(:,i) = X(:,i)/sqrt(norm(X(:,i),'fro')^2/M);
end
R_NS = (NX * NX'/N);
end

