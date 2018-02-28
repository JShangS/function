%%%%%%%%%%先感知rho的范围再估计协方差
clc
clear
close all
rou0=0.3;
rou=0.1:0.05:0.99;
str='p';
N=8;
L=round(2*N);
R0=fun_rho(rou0,N);
for i =1:length(rou)
    R=fun_rho(rou(i),N);
    MAM_t(:,:,i)=R;
end
%%%%%%%%%%%感知%%%%%%%%%
% X = fun_TrainData(str,N,L,R0,1,1,2);
% RX=fun_NSCMN(X);
% for i = 1:length(rou)
%     distance(i) = fun_ReimanDistance(RX,MAM_t(:,:,i));
% end
% rou_min = rou(find(min(distance)==distance));
% %%%%%%感知后的范围
% rou = max(0.1,rou_min-0.05):0.05:min(rou_min+0.05,0.99);
% %%%%%%%%%%%估计协方差
% for i =1:length(rou)
%     R=fun_rho(rou(i),N);
%     MAM(:,:,i)=R;
% end
LL=1000;
parfor i =1:LL
    x0 = fun_TrainData(str,N,1,R0,1,1,2);
    X = fun_TrainData(str,N,L,R0,1,1,2);
    RX=fun_NSCMN(X);
    Rx0=fun_SCM(x0);
    MAM = fun_congnition(X);
    [R_MAM,~,ratio]=fun_information_estimation(Rx0,MAM);
    error_R_MAM(i) = norm(R_MAM-R0,'fro')/norm(R0,'fro');%fun_ReimanDistance(R0,R_MAM);
    error_RX(i) = norm(RX-R0,'fro')/norm(R0,'fro');%fun_ReimanDistance(R0,RX);
end
m_errorR_MAM=mean(error_R_MAM);
m_errorRx0=mean(error_RX);
