function [ R_S ] = fun_SCM( X )
%����Э�������
%һ����һ�����뵥Ԫ
%X:ѵ������
[~,N] = size(X);
R_S = (X*X');
end

