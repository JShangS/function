function [ R_MAM, distance,ratio] = fun_information_estimation(R0, MAM)
%FUN_INFORMATION_ESTIMATION �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%%%%%������Ϣ���ζ�����Э������Ʒ�����
%%%R0:�ο�Э����
%%%MAM:��ģ�����Э����ΪN*N��LΪģ�͸���
[N,~,L]=size(MAM); 
distance = zeros(L,1);
for i=1:L
    distance(i) = fun_ReimanDistance(R0,MAM(:,:,i));
end
distance = distance/sum(distance);
Sum_ratio = sum(exp(-distance));
ratio = exp(-distance)/Sum_ratio;
ratio = reshape(ratio,1,1,L);
ratio = repmat(ratio,N,N,1);
R_MAM = sum(ratio.*MAM,3);
% ratio = ratio(1,1,1:L);
end

