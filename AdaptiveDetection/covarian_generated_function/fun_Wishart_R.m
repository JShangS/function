function [ RW ] = fun_Wishart_R(R,mu)
%FUN_WISHART_R �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%%���Ӿ�ֵΪR��wishart�ֲ���Э�����R����fun_rho�õ���Ҳ�����������Э���
%%<Conjugate Bayesian analysis of the Gaussian distribution>
%%%freedom���ɶ�
RW = wishrnd(1/mu*R,mu);
end

