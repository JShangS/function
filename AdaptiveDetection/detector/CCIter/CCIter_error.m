%%IterCCÿһ�εĵ������
clc
clear 
close all
% warning off
n = 2; %����������
str_train = 'p';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
sigma_t = 0.1;
opt = 'k';
rou = [0.10:0.10:0.90,0.95,0.99];  %%Э����������ɵĳ�������
L_rou = length(rou);
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
SNRout=-5:1:20; % ���SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% ϵͳ����ʸ��

iter = 100;
fornum = 1000;

error_iter = zeros(fornum,iter,L_rou);
alpha0 = zeros(fornum,iter,L_rou);
for i_rou = 1:L_rou
    i_rou
    R = fun_rho(rou(i_rou),N,2);
    R_KA = zeros(size(R));
    for i = 1:1000
        t = normrnd(1,sigma_t,N,1);%%0~0.5%%ʧ������
        R_KA = R_KA + R.*(t*t')/1000;
    end
    for i =1:fornum
    Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
    R_SCM = (fun_SCMN(Train));
    R_NSCM = (fun_NSCMN(Train));
    [R_CC,alpha0(i,1,i_rou)] = fun_CC(Train,R_SCM,R_KA);
    error_iter(i,1,i_rou) = norm(R_CC-R,'fro');
        for k = 1:iter-1
            [R_CC,alpha0(i,k+1,i_rou)] = fun_CC(Train,R_CC,R_KA);
            error_iter(i,k+1,i_rou) = norm(R_CC-R,'fro');
        end
    end
end
for i = 1:i_rou
    m_error_iter(i,:) = mean(error_iter(:,:,i),1)/norm(R,'fro');
    m_alpha0(i,:) = mean(alpha0(:,:,i),1);
end
figure(1)
hold on
str = [];
for i = 1:L_rou
plot([1:iter],m_error_iter(i,:));
str{i} = ['\rho=',num2str(rou(i))];
end
h_leg = legend(str);
xlabel('Number of iterations ','FontSize',10)
ylabel('NFN','FontSize',10)
set(gca,'FontSize',10)
set(h_leg,'Location','SouthEast')

figure(2)
hold on
for i = 1:L_rou
plot([1:iter],m_alpha0(i,:));
end
h_leg = legend(str);
xlabel('Number of iterations ','FontSize',10)
ylabel('\alpha','FontSize',10)
set(gca,'FontSize',10)
set(h_leg,'Location','SouthEast')

