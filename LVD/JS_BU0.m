function [ X ] = JS_BU0( x,offset,opt)
%JS_BU0 此处显示有关此函数的摘要
%   此处显示详细说明
%左移右移补0,上移下移补0函数
% x：向量或者数组
% offset ：移动量，正数是左移或上移，负数是右移下移
% opt: 1,左右移位，0上下移位，默认左右
if (nargin == 0)
    error('至少输入一个参数')
end
if (nargin == 1)
    offset = 0;
    opt = 1;
    disp('没做任何操作');
end
if (nargin == 2)
    opt = 1;
%     disp('默认左右移动');
end
if opt == 1  %左右
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
elseif opt == 0 %上下
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

