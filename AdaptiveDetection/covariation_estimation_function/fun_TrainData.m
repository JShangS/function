function [ X ] = fun_TrainData(N,L, M,opt)
%FUN_TRAINDATA �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%����ѵ������
%%N������ʸ��ά��
%%L��ѵ�����ݳ���
%%opt��������ѵ���������ͣ�1����˹��2��δ�����....
if nargin==3
    opt = 1;
end
if nargin<3
    error('�����������Ϊ3��')
end
%%������˹�Ӳ�
if opt == 1
    X = (randn(N,L)+1i*randn(N,L))/sqrt(2);
    M_half = M^0.5;
    X = M_half*X;
end
end

