%%%�����������Ȼ����뵽CC��,��AIWCM��������
clc
clear 
close all
Na = 4;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
rou = 0.95;  %%Э����������ɵĳ�������
rouR = zeros(N,N);  %%��ʵ���Ӳ�Э����
L=round(1.5*N); 
% theta_sig = 0.1;
% nn = 0:N-1;
% s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% ϵͳ����ʸ��
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);
    end
end
count = 1;
MC = 1000;
for i_MC =1:MC
    t = normrnd(1,0.2,N,1);%%0~0.5%%ʧ������
    R_KA = rouR.*(t*t');
    irouR=inv(rouR);
    rouR_half=rouR^0.5;
    ALL = sum(sum(sqrt(abs(rouR).^2)));%����Э����ķ���

    Train = fun_TrainData(N,L,rouR);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    x0 = fun_TrainData(N,1,rouR); % �����źŽ������Ӳ�������
    Sx0 = (x0*x0');

    beta_child = ones(L,1);
    beta_parent = zeros(L,1);
    R_MKA = R_KA;
    while (sum(abs(beta_parent-beta_child))/sum(abs(beta_child))>0.01)
        count;
        beta_parent = beta_child;
        R_KA_M_inv = inv(R_MKA);
        beta_child_t1 = Train'*R_KA_M_inv*x0;
        beta_child_t2 = diag(Train'*R_KA_M_inv*Train);
        beta_child = abs(beta_child_t1./beta_child_t2);
        R_KA_M_t = 0;
        for i = 1:L
            R_KA_M_t = R_KA_M_t+beta_child(i)*(Train(:,i)*Train(:,i)')/sum(beta_child);
        end
        R_MKA = abs(R_KA_M_t);
        count =count+1;
        if count >50
            break;
        end
     end
    error_MKA(i_MC) = sum(sum(sqrt(abs(R_MKA-rouR).^2)))/ALL;
    error_KA(i_MC) = sum(sum(sqrt(abs(R_KA-rouR).^2)))/ALL;
end
mean(error_MKA)
mean(error_KA)
var(error_MKA)
var(error_KA)

