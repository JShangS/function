function [ R_S ] = fun_SCMN( X )
%����Э�������/N
%һ����һ�����뵥Ԫ
%X:ѵ������
[~,N] = size(X);
R_S = (X*X'/N);
end

