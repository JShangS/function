function [ Tclglrt ] = fun_iGLCGLRT( lambda,mu,RKA,R,x0,p )
%FUN_CLGRT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%ɫ���ص�GLRT�Ľ��
[iR1,~] = fun_CLGLRT_icovariance(RKA,R,x0,p,1);%%H1�µ�Э������ƽ��
[iR0,~] = fun_CLGLRT_icovariance(RKA,R,x0,p,0);%%H0�µ�Э������ƽ��
[N,~]=size(x0);
a = (p'*iR1*x0)/(p'*iR1*p);
tmp1 = det(iR1)*((x0-a*p)'*iR1*(x0-a*p)+1/mu)^(-(N));%lambda
tmp0 = det(iR0)*(x0'*iR0*x0+1/mu)^(-(N));%lambda
Tclglrt = abs(tmp1/tmp0);
end

