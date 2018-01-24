function [ Tclglrt ] = fun_CLGLRT3(lambda,mu,RKA,R,x0,p )
%FUN_CLGRT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%ɫ���ص�GLRT�Ľ��
[R1,~] = fun_CLGLRT_covariance2(RKA,R,x0,p,1);%%H1�µ�Э������ƽ��
[R0,~] = fun_CLGLRT_covariance2(RKA,R,x0,p,0);%%H0�µ�Э������ƽ��
[N,~]=size(x0);
iR1 = inv(R1);
iR0 = inv(R0);
a = (p'*iR1*x0)/(p'*iR1*p);
tmp1 = det(R0)*((x0-a*p)'*iR1*(x0-a*p)+1/mu)^(-(N));%lambda
tmp0 = det(R1)*(x0'*iR0*x0+1/mu)^(-(N));%lambda
Tclglrt = abs(tmp1/tmp0);
end

