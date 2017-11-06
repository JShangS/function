function [ R_CC ] = fun_CC( X,R,R_KA )
%FUN_CC �˴���ʾ�йش˺�����ժҪ
%��On Using a priori Knowledge in Space-Time Adaptive Processing��
%   �˴���ʾ��ϸ˵��
%%ѵ���������Ƶ�Э���������Э�����������ϣ�����͹�Ż��õ����ϵ����
%%X:ѵ������
%R,�������Ƶ�Э����
%R_KA:����Э����
[M,N]=size(X);
rou_ba_t = 0;
% for i = 1:N
%     rou_ba_t = rou_ba_t+sum(abs(X(:,i))).^4/N^2;
% end
% rou_ba = rou_ba_t-sum(sum(abs(R).^2))/N;
rou_ba = sum(diag(X'*X).^2)/N^2-sum(sum(abs(R).^2))/N;%��18��ʽ,
alpha0 = rou_ba/(rou_ba+sum(sum(abs(R-R_KA).^2)));
R_CC = (1-alpha0)*R+alpha0*R_KA;
end

