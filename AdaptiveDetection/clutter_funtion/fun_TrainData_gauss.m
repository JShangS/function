function [ X ] = fun_TrainData_gauss(N,L,M)
%FUN_TRAINDATA �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%JerryShang��2017.11.09
%%����ѵ������
%%N������ʸ��ά��
%%L��ѵ�����ݳ���
if nargin<3
    error('�����������Ϊ3��')
end
%%������˹�Ӳ�
X = (randn(N,L)+1i*randn(N,L))/sqrt(2);
M_half = M^0.5;
X = M_half*X;
end

