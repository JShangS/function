function [ X ] = JS_BU0( x,offset,opt)
%JS_BU0 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%�������Ʋ�0,�������Ʋ�0����
% x��������������
% offset ���ƶ��������������ƻ����ƣ���������������
% opt: 1,������λ��0������λ��Ĭ������
if (nargin == 0)
    error('��������һ������')
end
if (nargin == 1)
    offset = 0;
    opt = 1;
    disp('û���κβ���');
end
if (nargin == 2)
    opt = 1;
%     disp('Ĭ�������ƶ�');
end
if opt == 1  %����
    [M,N] = size(x);
    if  abs(offset) <= N 
        if offset>=0
            X(:,1:N) = [x(:,1+offset:N),zeros(M,offset)];
        else
            X(:,1:N) = [zeros(M,abs(offset)),x(:,1:N-abs(offset))];
        end
    else
        X = zeros(M,N);
    end
elseif opt == 0 %����
    [M,N] = size(x);
    if  abs(offset) <= M 
        if offset>=0
            X(1:M,:) = [x(1+offset:M,:);zeros(abs(offset),N)];
        else
            X(1:M,:) = [zeros(abs(offset),N);x(1:M-abs(offset),:)];
        end
    else
        X = zeros(M,N);
    end
end
end

