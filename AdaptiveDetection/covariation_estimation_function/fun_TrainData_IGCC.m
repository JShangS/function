function [ Train ] = fun_TrainData_IGCC( N,L,M,lamda,mu,opt)
%FUN_TRAINDATA_IGCC �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%JerryShang��2017.11.16
%%��gamma�����ĸ��ϸ�˹�Ӳ�
%%N,L����,����ʸ��ά�ȣ��ο���Ԫ��Ŀ
%%M������˹Э����
%%lamda����״����
%%mu���߶Ȳ���
%%opt:1:ÿ�����뵥Ԫ������ֵ��һ����2��ÿ����Ԫ����ֵһ������ʱ�˻�ΪSIRP
if nargin ==5
    opt = 1;
end
if opt==1
    pct=sqrt(ones(L,1)./gamrnd(lamda,mu,L,1));
    X = fun_TrainData_gauss(N,L,M);
    Train=(ones(N,1)*pct.').*X;
elseif opt==2
    pct=sqrt(1./gamrnd(lamda,mu,1,1));
    X = fun_TrainData_gauss(N,L,M);
    Train=(ones(N,L)*pct.').*X;
else
    error('ֻ��1��2����ѡ��')
end

end
