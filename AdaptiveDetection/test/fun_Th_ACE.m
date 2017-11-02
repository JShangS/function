function [ Th_ACE ] = fun_Th_ACE( K,N,PFA )
%FUN_TH_ACE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%K��ѵ����������
%N����ⵥԪ����
%PFA:�龯��
t1 = 1-(PFA)^(1/(K-N+1));
t2 = 1-(K-N+1)/(K+1)*(PFA)^(1/(K-N+1));
Th_ACE = t1/t2;
end

