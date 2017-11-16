%  复合高斯
%  2014.06.11

clc;clear all;close all
N = 12;                     % 阵元数
K_p = 4;                    % 待测距离单元数
K_s=24;                   % 训练样本数
SNR_dB = [-10:2.5:20];    % 输入SNR [0:1.5:30]; 
SNR = 10.^( SNR_dB / 10 );  %  SNR数值


Mbar = zeros(N);              % 先验协方差矩阵
for ii = 1:N
    for jj = 1:N
    Mbar(ii,jj) = 0.9^(abs(ii-jj));  
    end
end
M=Mbar;  %  真实协方差

%=========================子空间参数=====================================%
r=1;
H=1/sqrt(N)*ones(N,r);  % N*r  信号导向矢量 
n=[0:N-1]';

%=========================目标信号幅度=====================================%
p_k_abs = (SNR/(H'*inv(Mbar)*H)/1); 
p_k_abs =sqrt(p_k_abs);  
%==============================================================%
Pfa=10^(-2);
runs=10^3;   % 统计次数 ：100/Pfa;
 
stas0_GLRT_ML=zeros(1,runs); 
stas1_GLRT_ML = zeros(runs,length(SNR_dB)); 
 

for oo=1:runs  % 完整的一次独立实验


pct=sqrt(ones(K_s,1)./gamrnd(8,1/5,K_s,1)); 
pct0=sqrt(ones(K_p,1)./gamrnd(8,1/5,K_p,1));
%------------------------------------------------------------------------
for ell=1:K_p
        y0(:,ell)=sqrtm(M)*(randn(N,1)+j*randn(N,1))/sqrt(2);
end
for e2=1:K_s
    y(:,e2)=sqrtm(M)*(randn(N,1)+j*randn(N,1))/sqrt(2); 
end     % 产生K_s个辅数据及采样协方差

     y=(ones(N,1)*pct.').*y;
     S = y*y';      % K_s采样协方差
     % 产生K_p个主数据
     y0=(ones(N,1)*pct0.').*y0;

%=========================ML estimate of M =====================================%
sum_0=0;
for j3=1:K_s
    sum_0=sum_0+y(:,j3)*(y(:,j3))'/((y(:,j3))'*y(:,j3));
end
M_ml0=(N/K_s)*sum_0;
for i3=1:3
    sum_2=0;
    iM_ml0=inv(M_ml0);
    for j4=1:K_s
        part3=y(:,j4)'*iM_ml0*y(:,j4);
        sum_2=sum_2+(y(:,j4)*(y(:,j4))')/part3;
    end
    M_ml0=(N/K_s)*sum_2;
end 
M_ml=M_ml0;
M_mmse=M_ml;
%=========================计算H0统计量======================================%
        iM_mmse=inv(M_mmse);  
        iM_ml=inv(M_ml);  
        Xp=y0;
                                 
          sum0_ml=0;
          for k1=1:K_p
              sum0_ml=sum0_ml+log(1-abs(H'*iM_ml*Xp(:,k1))^2/((H'*iM_ml*H)*(Xp(:,k1)'*iM_ml*Xp(:,k1)))); 
          end
         stas0_GLRT_ML(1,oo)=-N*sum0_ml;
%=========================计算H1统计量======================================%
   Xp1=zeros(N,K_p);
   for aa = 1:length(SNR )
       for bb=1:K_p
           phi0=2*pi*rand(1);
           p_k=p_k_abs*exp(1j*phi0); % 信号幅度 1*length(SNR)
           Xp1(:,bb)= y0(:,bb)+H*p_k(aa);
       end         
         sum1_ml=0;
          for k1=1:K_p
              sum1_ml=sum1_ml+log(1-abs(H'*iM_ml*Xp1(:,k1))^2/((H'*iM_ml*H)*(Xp1(:,k1)'*iM_ml*Xp1(:,k1)))); 
          end
          stas1_GLRT_ML(oo,aa)=-N*sum1_ml;
   end
end

%=========================计算阈值=====================================%

     stas0_GLRT_ML=sort(stas0_GLRT_ML,'descend');     % 统计量降序

     Th_GLRT_ML=stas0_GLRT_ML(10);
%=========================计算Pd=========================================%
 PD1= zeros(1,length(SNR_dB)); 

for ii=1:runs
    % H1一次统计量
    for aa = 1:length(p_k) 
        if stas1_GLRT_ML(ii,aa)>Th_GLRT_ML 
          PD1(1,aa)=PD1(1,aa)+1;              
        end  
    end
end
       
  PD1=PD1/runs;

%=========================绘图=========================================%
 plot(SNR_dB,PD1,'b--^','LineWidth',1.25);
 xlabel('SNR (dB)'); ylabel('Pd'); 
 legend('GLRT-ML');
grid on
% axis([5 35 0 1]);
% axis([0 30 0 1])