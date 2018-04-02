function [ X0 ] = fun_iterMAMH1( MAM,Z)
%FUN_ITERMAMH0 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%%��Adaptive detection of distributed targets in compound-Gaussian clutter 
%%without secondary data: An approach based on multiple a-priori spectral
%%models����H1�����µĹ���
%%MAM:��ģ��
%%Z����ⵥԪ����
iter_num = 50;%%��������
[N,L] = size(Z);
[~,~,num] = size(MAM);
X0 = zeros(N,N); %%%%��ʼЭ����
for i = 1:num
    X0 = X0 + 1/num * MAM(:,:,i);%%% ��ʼЭ����Ϊ��ģ�͵ľ�ֵ
end
tao = zeros(L,1);
T = zeros(N,N,L);
t = zeros(num,1);
for iter = 1:iter_num
    X0_t = zeros(N,N);
    for i = 1:L
        tao(i) = Z(:,i)' * Z(:,i) / L;
        for j = 1:L
             T(:,:,i) = T(:,:,i) + Z(:,j) * Z(:,j)' / tao(i);
        end
        t(i) = N / trace(X0 * T(:,:,i) / L);
        X0_t = X0_t + t(i) * MAM(:,:,i);
    end 
    X0 = X0_t / L;
end
end

