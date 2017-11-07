% 《Adaptively iterative weighting covariance matrix estimation 
% for airborne radar clutter suppression》，2015
%%%无先验知识，对已有的训练数据不在是平均加权而是根据和检测单元数据的关系进行加权

clc
clear 
close all
Na = 4;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
rou = 0.95;  %%协方差矩阵生成系数
rouR = zeros(N,N);
L=round(2*N); 
% theta_sig = 0.1;
% nn = 0:N-1;
% s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);
    end
end
t = normrnd(1,0.5,N,1);%%0~0.5%%失配向量
R_KA = rouR.*(t*t');
% rouR = R_KA;%diag(ones(1,N));%%%独立的协方差
irouR=inv(rouR);
rouR_half=rouR^0.5;
%%不完全Gamma函数
% % for i = 1:N
% %     for j = 1:L
% %         gamma_rand(i,j) = fun_IG(1,2); 
% %     end
% % end
X= (randn(N,L)+1i*randn(N,L))/sqrt(2);  % 产生方差为1的复高斯白噪声 % Rwhite1=1/snapshot1*X1*X1'; eig(Rwhite1); % round(mean(abs(eig(Rwhite1)))) == 1
% X= (gamma_rand.*randn(N,L)+1i*gamma_rand.*randn(N,L))/sqrt(2);
Train = rouR_half*X;%%产生的训练数据,协方差矩阵为rouR的高斯杂波
S=(rouR_half*X)*(rouR_half*X)'/L; % 有L个训练样本估计的杂波与噪声的协方差矩阵(Rhalf*X表示接收的L个训练数据)
S_abs = abs(S);
iS=inv(S);
W=(randn(N,1)+1i*randn(N,1))/sqrt(2); % 1i == -i
x0=rouR_half*W;%+pp; % noise=(randn(N,1)+j*randn(N,1))/sqrt(2);  % 接收信号仅包括杂波和噪声
R_x0 = abs(x0*x0');
Train_real = [Train]; %%文献的公式（35）
[~,Len] = size(Train_real);
beta_child = ones(Len,1);
beta_parent = zeros(Len,1);
R0_gu = eye(N,N);%R_x0;%eye(N,N);
count = 1;
while (norm(beta_parent-beta_child,'fro')/norm(beta_child,'fro')>0.1)
    count
    beta_parent = beta_child;
    R0_gu_inv = inv(R0_gu);
    beta_child_t1 = Train_real'*R0_gu_inv*x0;
    beta_child_t2 = diag(Train_real'*R0_gu_inv*Train_real);
    beta_child = abs(beta_child_t1./beta_child_t2).^2;
    beta_child_all(:,count) = beta_child;
    R0_gu_t = 0;
    for i = 1:Len
%         beta_child_t1 = Train(:,i)'*R0_gu_inv*x0;
%         beta_child_t2 = (Train(:,i)'*R0_gu_inv*Train(:,i));
%         beta_child(i) = (beta_child_t1/beta_child_t2).^2;
        R0_gu_t = R0_gu_t+beta_child(i)*(Train_real(:,i)*Train_real(:,i)')/sum(beta_child);%;sum(beta_child);
    end
    R0_gu = abs(R0_gu_t);
    count =count+1;
    if count >50
        break;
    end
end
erroS = norm(rouR - S_abs,'fro')/norm(rouR,'fro');
erro_AIWCM = norm(rouR - R0_gu,'fro')/norm(rouR,'fro');