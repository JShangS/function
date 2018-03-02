%%%协方差估计的信息几何方法，采用MAM模型为协方差的先验协方差。
clc
clear 
close all
rou0=0.95;
rou=0.9:0.01:0.99;
N=8;
L=round(1*N);
LL=1000;
R0=fun_rho(rou0,N);
for i =1:length(rou)
    R=fun_rho(rou(i),N);
    MAM(:,:,i)=R;
end

parfor i =1:LL
    x0 = fun_TrainData('k',N,1,R0,1,1,2);
    X = fun_TrainData('k',N,L,R0,1,1,2);
    RX=fun_NSCMN(X);
    Rx0=fun_SCM(x0);
    [R_MAM,~,ratio]=fun_information_estimation(RX,MAM);
    error_R_MAM(i) = norm(R_MAM-R0,'fro')/norm(R0,'fro');%fun_ReimanDistance(R0,R_MAM);
    error_RX(i) = norm(RX-R0,'fro')/norm(R0,'fro');%fun_ReimanDistance(R0,RX);.
end
m_errorR_MAM=mean(error_R_MAM);
m_errorRx0=mean(error_RX);

% distance1x0=0;
% distance2x0=0;
% for i =1:LL
%     i
%     X = fun_TrainData('k',N,L,R0,1,3,2);
%     R_x0=fun_SCMN(X);
% %     R_x0=fun_NSCMN(X);
%     distance1x0(i) = fun_ReimanDistance(R1,R_x0);
%     distance2x0(i) = fun_ReimanDistance(R2,R_x0);
%     ss=exp(-distance1x0(i))+exp(-distance2x0(i));
%     ratio1 = exp(-distance1x0(i))/ss;
%     ratio2 = exp(-distance2x0(i))/ss;
%     R_re=ratio1*R0+ratio2*R1;
%     error_R_re(i) = norm(R_re-R0,'fro')/norm(R0,'fro');
%     error_R_x0(i) = norm(R_x0-R0,'fro')/norm(R0,'fro');
% end
% m_R_re=mean(error_R_re);
% m_R_x0=mean(error_R_x0);
% m_distance1x0=mean(distance1x0);
% m_distance2x0=mean(distance2x0);

% x0 = fun_TrainData('p',N,1,R0,1,1,2);
% X = fun_TrainData('p',N,L,R0,1,1,2);
% RX=fun_NSCMN(X);
% Rx0=fun_NSCM(x0,1);
% for i = 1:length(rou)
%     distance1(i) = fun_ReimanDistance(Rx0,MAM(:,:,i));
%     distance2(i) = fun_ReimanDistance(RX,MAM(:,:,i));
% end
% plot(rou,distance1);
% figure
% plot(rou,distance2,'r');
% rou(find(min(distance2)==distance2))

