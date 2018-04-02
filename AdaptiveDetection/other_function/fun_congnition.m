function [ MAM,distance] = fun_congnition( X )
%FUN_CONGNITION �˴���ʾ�йش˺�����ժҪ
%%%%%%%%%%�ȸ�֪rho�ķ�Χ�ٹ���Э����
rou=0.1:0.05:0.99;
L_rou = length(rou);
[N,M] = size(X);
MAM_t = zeros(N,N,L_rou);
for i =1 : L_rou
    R=fun_rho(rou(i),N);
    MAM_t(:,:,i)=R;
end
if M<=N
    RX = fun_Positive(X,2);
else
    RX=fun_NSCMN(X);
end

distance = zeros(L_rou,1);
for i = 1 : L_rou
    distance(i) = fun_ReimanDistance(RX,MAM_t(:,:,i));
end
rou_min = rou(find(min(distance)==distance));
%%%%%%��֪��ķ�Χ
rou = max(0.1,rou_min-0.05):0.05:min(rou_min+0.05,0.99);
%%%%%%%%%%%����Э����
L_rou = length(rou);
MAM = zeros(N,N,L_rou);
for i =1:length(rou)
    R=fun_rho(rou(i),N);
    MAM(:,:,i)=R;
end
end

