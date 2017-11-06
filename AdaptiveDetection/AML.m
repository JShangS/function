% ��Covariance matrixestimation for CFAR detection in correlated heavy tailed clutter��
%%%������֪ʶ,��Ϊ��˹�ֲ�ʱ�������SCMһ��,
clc
clear 
close all
Na = 4;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
rou = 0.95;  %%Э�����������ϵ��
rouR = zeros(N,N);
L=round(4*N); 
% theta_sig = 0.1;
% nn = 0:N-1;
% s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% ϵͳ����ʸ��
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);
    end
end
irouR=inv(rouR);
rouR_half=rouR^0.5;
MC = 1000;
error_S = zeros(MC,1);
error_AML = zeros(MC,1);
error_NS = zeros(MC,1);
for iMC = 1:MC
    iMC
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
    while (abs(tao_child-tao_parent)>0.1)
%         count;
        tao_parent = tao_child;
        R0_gu_inv = inv(R0_gu);
        tao_child = diag(abs(Train_real'*R0_gu_inv*Train_real)/N);
    %     tao_child_all(count) = tao_child;
        R0_gu_t = 0;
        for i = 1:Len
            R0_gu_t = R0_gu_t+Train_real(:,i)*Train_real(:,i)'/Len/tao_child(i);
        end
        R0_gu = abs(R0_gu_t);
%         count =count+1;
%         if count >=3  
%             break;
%         end
    end
    error_S(iMC) = sum(sum(abs(rouR - S_abs)));
    error_AML(iMC) = sum(sum(abs(rouR - R0_gu)));
    error_NS(iMC) = sum(sum(abs(rouR - NS_abs)));
end
mean_error_S = mean(error_S);
mean_error_AML = mean(error_AML);
mean_error_NS = mean(error_NS);