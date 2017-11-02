% ��AML�Ľ����CC��͹�Ż����Ľ�����
%%1.��Covariance matrixestimation for CFAR detection in correlated heavy tailed clutter��
%%��2.<On Using a priori Knowledge in Space-Time Adaptive Processing>
%%%������֪ʶ
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
t = 1*rand(N,1);
R_KA = rouR.*(t*t'); %����Э����
irouR=inv(rouR);
rouR_half=rouR^0.5;
%%����ѵ������
X = (randn(N,L)+1i*randn(N,L))/sqrt(2);  % ��������Ϊ1�ĸ���˹������ % Rwhite1=1/snapshot1*X1*X1'; eig(Rwhite1); % round(mean(abs(eig(Rwhite1)))) == 1
Train = rouR_half*X;%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
%%SCM
S = (rouR_half*X)*(rouR_half*X)'/L; % ��L��ѵ���������Ƶ��Ӳ���������Э�������(Rhalf*X��ʾ���յ�L��ѵ������)
S_abs = abs(S);
iS=inv(S);
%%NSCM
for i = 1:L
    NX(:,i) = Train(:,i)/sqrt(Train(:,i)'*Train(:,i)/N);
end
NS = NX * NX'/L;
NS_abs = abs(NS);
%%��ⵥԪ
W=(randn(N,1)+1i*randn(N,1))/sqrt(2); % 1i == -i
x0=rouR_half*W;%+pp; % noise=(randn(N,1)+j*randn(N,1))/sqrt(2);  % �����źŽ������Ӳ�������
R_x0 = abs(x0*x0');
%%AML
Train_real = [Train]; 
[~,Len] = size(Train_real);
R0_gu = eye(N,N);%R_x0;%eye(N,N);%�Ե�λ��Ϊ������ֵ�ǣ��ڶ��ε������ΪNSCM���
tao_child = 1;%%���ε���ֵ
tao_parent = 0;%%�ϴε���ֵ
count = 1;
while (1)
    count
    tao_parent = tao_child;
    R0_gu_inv = inv(R0_gu);
    tao_child = abs(Train_real'*R0_gu_inv*Train_real)/N;
%     tao_child_all(count) = tao_child;
    R0_gu_t = 0;
    for i = 1:Len
        R0_gu_t = R0_gu_t+Train_real(:,i)*Train_real(:,i)'/Len;
    end
    R0_gu = (R0_gu_t);
    count =count+1;
    if count >=3  
        break;
    end
end
%%%AML+CC
%%CC
% sum_t = 0;
% for i=1:L
%     sum_t = sum_t + (Train_real(:,L)'*Train_real(:,L))^2;
% end
%%%AML+CC
rou_ba = sum(diag(Train_real'*Train_real).^2)/L^2-sum(sum(R0_gu.^2))/L;%2�ģ�18��ʽ,
alpha0 = rou_ba/(rou_ba+sum(sum((R0_gu-R_KA).^2)));
R_AML_CC = (1-alpha0)*R0_gu+alpha0*R_KA;

%%%S+CC
rou_ba_S = sum(diag(Train_real'*Train_real).^2)/L^2-sum(sum(S.^2))/L;%2�ģ�18��ʽ,
alpha0_S = rou_ba/(rou_ba+sum(sum((S-R_KA).^2)));
R_S_CC = (1-alpha0_S)*S+alpha0_S*R_KA;
%%%%
error_S = sum(sum(abs(rouR - S_abs)));
error_AML = sum(sum(abs(rouR - R0_gu)));
error_NS = sum(sum(abs(rouR - NS_abs)));
error_AML_CC = sum(sum(abs(rouR - R_AML_CC)));
error_S_CC = sum(sum(abs(rouR - R_S_CC)));