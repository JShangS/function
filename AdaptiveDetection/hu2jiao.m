function [ y ] = hu2jiao( x )
%HU2JIAO 此处显示有关此函数的摘要
%   弧度转角度
y = mod(x,pi)/pi*180;
end

