function [ R_CLGRT ] = fun_CLGLRT_covariance( lambda,mu,RKA,R,x0,p,opt )
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
   alpha1 = abs(lambda*mu*((x0-a1*p)'*iRKA*(x0-a1*p))/N);
   %%%beta1
   a2 = (p'*iR*x0)/(p'*iR*p);
   beta1 = abs(lambda*mu*((x0-a2*p)'*iR*(x0-a2*p))/N);
   %%%R_CLGRT
   sum_t1 = alpha1 + beta1;
   R_CLGRT = (alpha1/sum_t1)*RKA+(beta1/sum_t1)*R;
%     R_CLGRT = alpha1*RKA + beta1*R;
elseif opt == 0%H0�����µ�Э�������
    %%%alpha0
    alpha0 = abs(lambda*mu*(x0'*iRKA*x0)/N);
    %%%beta0
    beta0 = abs(lambda*mu*(x0'*iR*x0)/N);
     %%%R_CLGRT
    sum_t0 = alpha0 + beta0;
    R_CLGRT = (alpha0/sum_t0)*RKA+(beta0/sum_t0)*R;
%     R_CLGRT = alpha0*RKA + beta0*R;
end
end

