function [ MAM,distance] = fun_congnition( X )
%FUN_CONGNITION 此处显示有关此函数的摘要
%%%%%%%%%%先感知rho的范围再估计协方差
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
%%%%%%感知后的范围
rou = max(0.1,rou_min-0.05):0.05:min(rou_min+0.05,0.99);
%%%%%%%%%%%估计协方差
L_rou = length(rou);
MAM = zeros(N,N,L_rou);
for i =1:length(rou)
    R=fun_rho(rou(i),N);
    MAM(:,:,i)=R;
end
end

