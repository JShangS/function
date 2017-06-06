function [ f ] = f1( x1,x2 )
%F1 此处显示有关此函数的摘要
%   此处显示详细说明
%f=-(1+x.*sin(4*pi*x)-y.*sin(4*pi*y+pi)), -1<=x,y<=1
f = -(1+x1.*sin(4*pi*x1)-x2.*sin(4*pi*x2+pi));

end

