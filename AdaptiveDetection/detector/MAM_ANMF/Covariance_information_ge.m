%%%协方差估计的信息几何方法，采用MAM模型为协方差的先验协方差。
clc
clear 
close all
rou0=0.95;
% rou=0.5:0.05:0.99;
rou=[0.1,0.6,0.90];
N=8;
L=round(2*N);
LL=1000;
R0=fun_rho(rou0,N);
for i =1:length(rou)
    R=fun_rho(rou(i),N);
    MAM(:,:,i)=R;
end
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'/sqrt(N); %%%%%% 系统导向矢量
tic
for i =1:LL
    i
    x0 = fun_TrainData('p',N,1,R0,2,1,1);
    x0 = 10*s + x0;
    X = fun_TrainData('p',N,L,R0,2,1,1);
    RNX = abs(fun_NSCMN(X));
    RX = abs(fun_SCMN(X));
%     Rx0=fun_SCMN(x0);
    Rx0 = x0*x0';
    [V,D] = svd(Rx0);
    Evalue = abs(diag(D));
    index_1 = find(Evalue<=1);
    Evalue(index_1) = 1;
    % Evalue = sort(Evalue,'descend');
    D = diag(Evalue);
    Rx0 = abs(V)*D*abs(V');
    
%     [R_MAM_r,~,ratio_r]=fun_information_estimation(Rx0,MAM,'r');
    [R_MAM_c,~,ratio_c]=fun_information_estimation(Rx0,MAM,'c');
    [R_MAM_e,~,ratio_e]=fun_information_estimation(Rx0,MAM,'e');
    [R_MAM_l,~,ratio_l]=fun_information_estimation(Rx0,MAM,'l');
%     [R_MAM_p,~,ratio_p]=fun_information_estimation(Rx0,MAM,'p');
%     [R_MAM_ro,~,ratio_ro]=fun_information_estimation(Rx0,MAM,'ro');
    
%     error_R_MAM_r(i) = norm(R_MAM_r-R0,'fro')/norm(R0,'fro');%fun_ReimanDistance(R0,R_MAM);
    error_R_MAM_c(i) = norm(R_MAM_c-R0,'fro')/norm(R0,'fro');
    error_R_MAM_e(i) = norm(R_MAM_e-R0,'fro')/norm(R0,'fro');
    error_R_MAM_l(i) = norm(R_MAM_l-R0,'fro')/norm(R0,'fro');
%     error_R_MAM_p(i) = norm(R_MAM_p-R0,'fro')/norm(R0,'fro');
%     error_R_MAM_ro(i) = norm(R_MAM_ro-R0,'fro')/norm(R0,'fro');
    error_RNX(i) = norm(RNX-R0,'fro')/norm(R0,'fro');%fun_ReimanDistance(R0,RX);.
    error_RX(i) = norm(RX-R0,'fro')/norm(R0,'fro');%fun_ReimanDistance(R0,RX);.
end
toc
% m_errorR_MAM_r=mean(error_R_MAM_r)*100;
m_errorR_MAM_c=mean(error_R_MAM_c)*100;
m_errorR_MAM_e=mean(error_R_MAM_e)*100;
m_errorR_MAM_l=mean(error_R_MAM_l)*100;
% m_errorR_MAM_p=mean(error_R_MAM_p)*100;
% m_errorR_MAM_ro=mean(error_R_MAM_ro)*100;
m_errorRNX=mean(error_RNX)*100;
m_errorRX=mean(error_RX)*100;
str_bar = {'Euclidean','LogEuclidean','Cholesky','NSCM','SCM'};
% str_bar = {'Reiman','Euclidean','LogEuclidean',...
%     'Power-Euclidean','RootEuclidean','Cholesky','NSCM','SCM'};
data = [m_errorR_MAM_e,m_errorR_MAM_l,m_errorR_MAM_c,m_errorRNX,m_errorRX];
b = bar(data,'stack');
set(gca,'xticklabel',str_bar);

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

