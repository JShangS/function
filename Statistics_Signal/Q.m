function [ y ] = Q( x )
%   ��˹��β����
%   xΪ�������ޣ�yΪ����ֵ����׼��̬�ֲ�X>x�ĸ��� 
y = 0.5*erfc(x/sqrt(2));
end

