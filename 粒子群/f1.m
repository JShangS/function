function [ f ] = f1( x1,x2 )
%F1 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%f=-(1+x.*sin(4*pi*x)-y.*sin(4*pi*y+pi)), -1<=x,y<=1
f = -(1+x1.*sin(4*pi*x1)-x2.*sin(4*pi*x2+pi));

end

