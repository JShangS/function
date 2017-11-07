% ��Adaptively iterative weighting covariance matrix estimation 
% for airborne radar clutter suppression����2015
%%%������֪ʶ�������е�ѵ�����ݲ�����ƽ����Ȩ���Ǹ��ݺͼ�ⵥԪ���ݵĹ�ϵ���м�Ȩ

clc
clear 
close all
Na = 4;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
rou = 0.95;  %%Э�����������ϵ��
rouR = zeros(N,N);
L=round(2*N); 
% theta_sig = 0.1;
% nn = 0:N-1;
% s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% ϵͳ����ʸ��
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);
    end
end
t = normrnd(1,0.5,N,1);%%0~0.5%%ʧ������
R_KA = rouR.*(t*t');
% rouR = R_KA;%diag(ones(1,N));%%%������Э����
irouR=inv(rouR);
rouR_half=rouR^0.5;
%%����ȫGamma����
% % for i = 1:N
% %     for j = 1:L
% %         gamma_rand(i,j) = fun_IG(1,2); 
% %     end
% % end
X= (randn(N,L)+1i*randn(N,L))/sqrt(2);  % ��������Ϊ1�ĸ���˹������ % Rwhite1=1/snapshot1*X1*X1'; eig(Rwhite1); % round(mean(abs(eig(Rwhite1)))) == 1
% X= (gamma_rand.*randn(N,L)+1i*gamma_rand.*randn(N,L))/sqrt(2);
Train = rouR_half*X;%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
S=(rouR_half*X)*(rouR_half*X)'/L; % ��L��ѵ���������Ƶ��Ӳ���������Э�������(Rhalf*X��ʾ���յ�L��ѵ������)
S_abs = abs(S);
iS=inv(S);
W=(randn(N,1)+1i*randn(N,1))/sqrt(2); % 1i == -i
x0=rouR_half*W;%+pp; % noise=(randn(N,1)+j*randn(N,1))/sqrt(2);  % �����źŽ������Ӳ�������
R_x0 = abs(x0*x0');
Train_real = [Train]; %%���׵Ĺ�ʽ��35��
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