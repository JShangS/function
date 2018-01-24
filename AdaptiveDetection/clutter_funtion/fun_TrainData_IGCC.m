function [ Train ] = fun_TrainData_IGCC( N,L,M,lambda,mu,opt_train)
%FUN_TRAINDATA_IGCC �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%JerryShang��2017.11.16
%%��gamma����ĸ��ϸ�˹�Ӳ�
%%N,L����,����ʸ��ά�ȣ��ο���Ԫ��Ŀ
%%M������˹Э����
%%lamda����״����
%%mu���߶Ȳ���
%%opt:1: ÿ�����뵥Ԫ������ֵ��һ������CPI��ÿ�����뵥Ԫ������ֵһ����SIRP
%%%%%%2��ÿ����Ԫ����ֵһ����Ϊ���־��Ȼ���
%%%%%%3��ÿ����Ԫ����ֵÿʱÿ�̶���һ����������compound Gaussian
%%%�������������
if nargin == 5
    opt_train = 1;
end
if opt_train == 1
    pct=sqrt(ones(L,1)./gamrnd(lambda,mu,L,1));
    X = fun_TrainData_gauss(N,L,M);
    Train=(ones(N,1)*pct.').*X;
elseif opt_train == 2
    pct=sqrt(1./gamrnd(lambda,mu,1,1));
    X = fun_TrainData_gauss(N,L,M);
    Train=(ones(N,L)*pct.').*X;
elseif opt_train == 3
    pct = zeros(N,L);
    for i = 1:L
        pct(:,i)=sqrt(1./gamrnd(lambda,mu,N,1));
    end  
    X = fun_TrainData_gauss(N,L,M);
    Train=pct.*X;
else
    error('ֻ��1��2, 3����ѡ��')
end

end

