function [ Train ] = fun_TrainData_IGCC( N,L,M,lamda,mu )
%FUN_TRAINDATA_IGCC �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%��gamma����ĸ��ϸ�˹�Ӳ�
%%N,L����,����ʸ��ά�ȣ��ο���Ԫ��Ŀ
%%M������˹Э����
%%lamda����״����
%%mu���߶Ȳ���
pct=sqrt(ones(L,1)./gamrnd(lamda,mu,L,1));
X = fun_TrainData_gauss(N,L,M);
Train=(ones(N,1)*pct.').*X;
end

