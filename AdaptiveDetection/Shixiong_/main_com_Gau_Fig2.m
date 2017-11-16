clear ;
close all;
clc;
%《Geometric barycenters for covariance estimation in compound-Gaussian clutter》
 
N=8;   %% the number of pulse in a CPI
K=10;   %% the number of secondary data

fd=0.15; %% the doppler of target
P=exp(1j*2*pi*fd*(0:1:N-1)).'; %% streer vector

%% =====Covariance matrix of clutter====%%%%%
v=0.5;  %% the parameters of the clutter
rou=0.99;  % one-lag correlation coefficient 
u=1;

fdc=0.05;  % clutter normalised doppler freqquency
deltdB=20; % CNR=20dB
delt=10.^(deltdB/10); 
Sigma0=zeros(N,N); % covariance matrix of clutter
for ii=1:N
    for k=1:N
         Sigma0(ii,k)=delt*(rou^(abs(ii-k)))*exp(1j*2*pi*fdc*(ii-k));
    end
end

Sigma=zeros(N,N); % matrix of disturbance
Sigma=Sigma0+eye(N); % clutter plus noise

[E,D]=eig(Sigma);    % Transform matrix 
RT=E*D.^0.5;%%%

eig_valu=diag(D);
eig_valu_max=max(eig_valu);
eig_valu_min=min(eig_valu);
K_m=eig_valu_max/eig_valu_min;

%% =====threshold calculation=====%%%%
pfa=10^(-4);
num_pfa=100/pfa; % Number of Monte_Carlo simulation
vet_GLRT=zeros(1,num_pfa); 
vet_NSCM=zeros(1,num_pfa); 
vet_FPE=zeros(1,num_pfa);
vet_E=zeros(1,num_pfa);  
vet_L=zeros(1,num_pfa);
vet_H=zeros(1,num_pfa);
vet_A=zeros(1,num_pfa);
vet_C=zeros(1,num_pfa);

for h=1:num_pfa
     h
      %%%%=== generate the secondary data========%%%%%
      %%%%纹理是个K分布的纹理 
      Snd_r=zeros(N,K);
      for k=1:K
            gs=sqrt(1/2)*(randn(N,1)+1j*randn(N,1)); 
            taos=gamrnd(v,1/v,1,1);  %K-distribution clutter
            Snd_r(:,k)=sqrt(taos)*(RT*gs);
      end
      
      %%===== Estimation of covariance matrix with traditional  methods===%%  
      Re_NSCM=NSCM(Snd_r);  %% NSCM
      Re_FPE=FPE(Snd_r,4);     %% FPE
      
      %%===== Estimation of covariance matrix with geometric barycenters methods====%% 
      Sk=zeros(N,N,K);
      Sk=sample_matrix_Norma(Snd_r,K_m); %% produce Sk
      
      Re_E=stima_euclidea(Sk);                     %% Euclidean estimator
      Re_L=stima_log_euclidea(Sk);               %% Log-Euclidean estimator
      Re_H=stima_root_euclidea(Sk);             %% Root-Euclidean estimator
      Re_A=stima_power_euclidea(Sk,0.8);    %% Power-Euclidean estimator
      Re_C=stima_cholesky(Sk);                    %% Cholesky estimator
            
      %%%%===generate the received data of the CUT=======%%%%%%%%%
      r=zeros(N,1);  
      g=sqrt(1/2)*(randn(N,1)+j*randn(N,1)); %k-distribution clutter
      tao=gamrnd(v,1/v,1,1);
      r=sqrt(tao)*(RT*g);
      
      %%%====test design
      vet_GLRT(1,h)=((abs(P'*inv(Sigma)*r))^2)/((P'*inv(Sigma)*P)*(r'*inv(Sigma)*r));
      vet_NSCM(1,h)=((abs(P'*inv(Re_NSCM)*r))^2)/((P'*inv(Re_NSCM)*P)*(r'*inv(Re_NSCM)*r));
      vet_FPE(1,h)=((abs(P'*inv(Re_FPE)*r))^2)/((P'*inv(Re_FPE)*P)*(r'*inv(Re_FPE)*r));
      vet_E(1,h)= ((abs(P'*inv(Re_E)*r))^2)/((P'*inv(Re_E)*P)*(r'*inv(Re_E)*r));
      vet_L(1,h)=((abs(P'*inv(Re_L)*r))^2)/((P'*inv(Re_L)*P)*(r'*inv(Re_L)*r));
      vet_H(1,h)=((abs(P'*inv(Re_H)*r))^2)/((P'*inv(Re_H)*P)*(r'*inv(Re_H)*r));
      vet_A(1,h)=((abs(P'*inv(Re_A)*r))^2)/((P'*inv(Re_A)*P)*(r'*inv(Re_A)*r));
      vet_C(1,h)=((abs(P'*inv(Re_C)*r))^2)/((P'*inv(Re_C)*P)*(r'*inv(Re_C)*r));
      
end

vet_GLRT=sort(vet_GLRT);
Threshold_GLRT=vet_GLRT(num_pfa-100); % GLRT Theshold
vet_NSCM=sort(vet_NSCM);
Threshold_NSCM=vet_NSCM(num_pfa-100);
vet_FPE=sort(vet_FPE);
Threshold_FPE=vet_FPE(num_pfa-100);
vet_E=sort(vet_E);
Threshold_E=vet_E(num_pfa-100);
vet_L=sort(vet_L);
Threshold_L=vet_L(num_pfa-100);
vet_H=sort(vet_H);
Threshold_H=vet_H(num_pfa-100);
vet_A=sort(vet_A);
Threshold_A=vet_A(num_pfa-100);
vet_C=sort(vet_C);
Threshold_C=vet_C(num_pfa-100);

%% %%%%==== Pd=====%%%%%
tic
SINR=linspace(-10,40,21);
num_pd=6000;
pd_GLRT=zeros(size(SINR));
pd_NSCM=zeros(size(SINR));
pd_FPE=zeros(size(SINR));
pd_E=zeros(size(SINR));
pd_L=zeros(size(SINR));
pd_H=zeros(size(SINR));
pd_A=zeros(size(SINR));
pd_C=zeros(size(SINR));

for ll1=1:num_pd
     ll1
    for ll2=1:length(SINR)
               
          %%%%=== generate the secondary data========%%%%%
          Snd_r=zeros(N,K);
          for k=1:K
                gs=sqrt(1/2)*(randn(N,1)+j*randn(N,1)); 
                taos=gamrnd(v,1/v,1,1);  %K-distribution clutter
                Snd_r(:,k)=sqrt(taos)*(RT*gs);
          end

          %%===== Estimation of covariance matrix with traditional  methods===%%  
          Re_NSCM=NSCM(Snd_r);  %% NSCM
          Re_FPE=FPE(Snd_r,4);     %% FPE

          %%===== Estimation of covariance matrix with geometric barycenters methods====%% 
          Sk=zeros(N,N,K);
          Sk=sample_matrix_Norma(Snd_r,K_m); %% produce Sk

          Re_E=stima_euclidea(Sk);                     %% Euclidean estimator
          Re_L=stima_log_euclidea(Sk);               %% Log-Euclidean estimator
          Re_H=stima_root_euclidea(Sk);             %% Root-Euclidean estimator
          Re_A=stima_power_euclidea(Sk,0.8);    %% Power-Euclidean estimator
          Re_C=stima_cholesky(Sk);     
                  
          %%%%===generate the received data of the CUT=======%%%%%%%%%
          %%%generate the clutter  of the CUT
          rc=zeros(N,1);  
          gc=sqrt(1/2)*(randn(N,1)+1j*randn(N,1)); %k-distribution clutter
          taoc=gamrnd(v,1/v,1,1);
          rc=sqrt(taoc)*(RT*gc);
          %%%SNIR of steer vector
          Tamplite=(10^(SINR(ll2)/10))/(P'*inv(Sigma)*P);
%           Talfai=sqrt(Tamplite/2)*(randn(1,1)+j*randn(1,1)); %%  target fluctuation
          r=sqrt(Tamplite)*exp(1j*2*pi*rand(1))*P+rc; %% recived signal
          
          %%%====test design
          T_GLRT=((abs(P'*inv(Sigma)*r))^2)/((P'*inv(Sigma)*P)*(r'*inv(Sigma)*r));
          T_NSCM=((abs(P'*inv(Re_NSCM)*r))^2)/((P'*inv(Re_NSCM)*P)*(r'*inv(Re_NSCM)*r));
          T_FPE=((abs(P'*inv(Re_FPE)*r))^2)/((P'*inv(Re_FPE)*P)*(r'*inv(Re_FPE)*r));
          T_E= ((abs(P'*inv(Re_E)*r))^2)/((P'*inv(Re_E)*P)*(r'*inv(Re_E)*r));
          T_L=((abs(P'*inv(Re_L)*r))^2)/((P'*inv(Re_L)*P)*(r'*inv(Re_L)*r));
          T_H=((abs(P'*inv(Re_H)*r))^2)/((P'*inv(Re_H)*P)*(r'*inv(Re_H)*r));
          T_A=((abs(P'*inv(Re_A)*r))^2)/((P'*inv(Re_A)*P)*(r'*inv(Re_A)*r));
          T_C=((abs(P'*inv(Re_C)*r))^2)/((P'*inv(Re_C)*P)*(r'*inv(Re_C)*r));
          
          %%Calculation of pd
       if T_GLRT>Threshold_GLRT
           pd_GLRT(ll2)=pd_GLRT(ll2)+1;
       end 
       if T_NSCM>Threshold_NSCM
           pd_NSCM(ll2)=pd_NSCM(ll2)+1;
       end 
       if T_FPE>Threshold_FPE
           pd_FPE(ll2)=pd_FPE(ll2)+1;
       end
       if T_E>Threshold_E
           pd_E(ll2)=pd_E(ll2)+1;
       end
        if T_L>Threshold_L
           pd_L(ll2)=pd_L(ll2)+1;
        end
        if T_H>Threshold_H
           pd_H(ll2)=pd_H(ll2)+1;
        end
        if T_A>Threshold_A
           pd_A(ll2)=pd_A(ll2)+1;
        end
         if T_C>Threshold_C
           pd_C(ll2)=pd_C(ll2)+1;
        end
        
    end
end
toc

save PDN8K10V05ROU099


labeltsize=12;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 7;
figure;axes1 = axes;
plot(SINR,pd_GLRT/num_pd,'k-','LineWidth',linewide1,'MarkerSize',mkft);hold on;
plot(SINR,pd_FPE/num_pd,'r-.+','LineWidth',linewide1,'MarkerSize',mkft);hold on;
plot(SINR,pd_NSCM/num_pd,'k-.o','LineWidth',linewide1,'MarkerSize',mkft);hold on;
plot(SINR,pd_E/num_pd,'r--*','LineWidth',linewide1,'MarkerSize',mkft);hold on;
plot(SINR,pd_L/num_pd,'b-^','LineWidth',linewide1,'MarkerSize',mkft);hold on;
plot(SINR,pd_H/num_pd,'g-o','LineWidth',linewide1,'MarkerSize',mkft);hold on;
plot(SINR,pd_A/num_pd,'m--s','LineWidth',linewide1,'MarkerSize',mkft);hold on;
plot(SINR,pd_C/num_pd,'c--+','LineWidth',linewide1,'MarkerSize',mkft);hold on;
title('N=8, K=10, \nu=0.5, \rho=0.99');
legend1 = legend(axes1,{'NMF','FPE','NSCM','Euclidean','Log-Euclidean','Root-Euclidean','Power-Euclidean', 'Cholesky'});
%set(legend1,'Location','South East');
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn);
xlabel('SINR','FontSize',labeltsize,'FontWeight',fw,'FontName',fn);