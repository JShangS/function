function [ Train ] = fun_TrainData_K( N,L, M, v)
%FUN_TRAINDATA_K �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%K����ֲ��ĸ��ϸ�˹�ֲ�
%%N,L����
%%M������˹Э����
[U,D] = eig(M);
M_half = U*D.^0.5;
Train=zeros(N,L);
for l=1:L
    gs=sqrt(1/2)*(randn(N,1)+1j*randn(N,1)); 
    taos=gamrnd(v,1/v,1,1);  %K-distribution clutter
    Train(:,l)=sqrt(taos)*(M_half*gs);
end
end

