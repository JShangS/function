function [ Th_AMF ] = fun_Th_AMF( K,N,PFA )
%FUN_TH_ACE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%K��ѵ����������
%N����ⵥԪ����
%PFA:�龯��
Th_AMF = (K+1)/(K-N+1)*((PFA)^(1/(K-N+2))-1)
end

