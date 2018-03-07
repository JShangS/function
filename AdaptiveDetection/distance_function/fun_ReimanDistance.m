function [ distance ] = fun_ReimanDistance( A,B )
%FUN_REIMANDISTANCE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%%����������������ξ��룬��������<������Ϣ���ε��״�Ŀ����_�ڷ���>��ʽ4.18��4.20���
T = A^(-0.5)*B*A^(-0.5);
[V,D] = eig(T);%%%�����ֽ⣬VΪ����������DΪ����ֵ�Խ���
lnD = abs(log(D).^2);
distance = trace(lnD);
end

