function [ R_CLGRT,alpha,beta ] = fun_CLGLRT_covariance2(RKA,R,x0,p,opt )
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
   a1 = (p'*iRKA*x0)/(p'*iRKA*p);
   a2 = (p'*iR*x0)/(p'*iR*p);
   fai1 = (((x0-a1*p)'*iRKA*(x0-a1*p))/((x0-a2*p)'*iR*(x0-a2*p)));%fai1 = (((x0-a1*p)'*iRKA*(x0-a1*p))/((x0-a2*p)'*iR*(x0-a2*p)));
   alpha1 = (1/(1+fai1));
   beta1 = (fai1/(1+fai1));
   %%%R_CLGRT
   R_CLGRT = alpha1*RKA+beta1*R;
   alpha = alpha1;
   beta = beta1;
%     R_CLGRT = alpha1*RKA + beta1*R;
elseif opt == 0%H0�����µ�Э�������
    fai0 =  ((x0'*iRKA*x0)/(x0'*iR*x0));%fai0 =  ((x0'*iRKA*x0)/(x0'*iR*x0));
    alpha0 = (1/(1+fai0));
    beta0 = (fai0/(1+fai0));
    R_CLGRT = alpha0*RKA+beta0*R;
    alpha = alpha0;
    beta = beta0;
end
end

