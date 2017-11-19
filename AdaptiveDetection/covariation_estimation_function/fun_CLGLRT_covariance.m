function [ R_CLGRT ] = fun_CLGLRT_covariance( lamda,mu,RKA,R,x0,p,opt )
%FUN_CLGLRT_CON �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%lamda:��gamma�ֲ�����״������mu���߶Ȳ���
%RKA������Э���R:����Э����
%x0������ⵥԪ
%opt��H1��H0�����µ�Э������ƽ��
%p������ʸ��
iRKA = inv(RKA);
iR = inv(R);
[N,~] = size(x0);
if opt == 1 %%H1�����µ�Э�������
   %%%alpha1
   a1 = (p'*iRKA*x0)/(p'*iRKA*p);
   tmp11 = (lamda+N)*iRKA;
   tmp12 = N*((x0-a1*p)'*iRKA*(x0-a1*p));
   alpha1 = mu*(tmp11-tmp12)/N;
   %%%beta1
   a2 = (p'*iR*x0)/(p'*iR*p);
   tmp11 = (lamda+N)*iR;
   tmp12 = N*((x0-a2*p)'*iR*(x0-a2*p));
   beta1 = mu*(tmp11-tmp12)/N;
   %%%R_CLGRT
   R_CLGRT = alpha1*RKA+beta1*R;
elseif opt == 0%H0�����µ�Э�������
    %%%alpha0
    tmp01 = (lamda+N)*iRKA;
    tmp02 = N*(x0'*iRKA*x0);
    alpha0 = mu*(tmp01-tmp02)/N;
    %%%beta0
    tmp01 = (lamda+N)*iR;
    tmp02 = N*(x0'*iR*x0);
    beta0 = mu*(tmp01-tmp02)/N;
    R_CLGRT = alpha0*RKA+beta0*R;
end
end

