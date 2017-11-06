%%%对先验的修正然后加入到CC中,用AML进行修正
clc
clear 
close all
Na = 4;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
rou = 0.95;  %%协方差矩阵生成的迟滞因子
rouR = zeros(N,N);  %%真实的杂波协方差
L=round(2*N); 
% theta_sig = 0.1;
% nn = 0:N-1;
% s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);
    end
end
count = 1;
MC = 1000;
for i_MC =1:MC
    i_MC
    t = normrnd(1,0.3,N,1);%%0~0.5%%失配向量
    R_KA = rouR.*(t*t');
    irouR=inv(rouR);
    rouR_half=rouR^0.5;
    ALL = sum(sum(sqrt(abs(rouR).^2)));%真正协方差的范数

    Train = fun_TrainData(N,L,rouR);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = fun_TrainData(N,1,rouR); % 接收信号仅包括杂波和噪声
    Sx0 = (x0*x0');

    tao_child = 1;
    tao_parent = 0;
    R_MKA = R_KA;
    while (abs(tao_child-tao_parent)>0.1)%
        tao_parent = tao_child;
        R_MKA_inv = inv(R_MKA);
        tao_child = diag(abs(Train'*R_MKA_inv*Train)/N);
        R_MKA_t = 0;
        for i = 1:L
            R_MKA_t = R_MKA_t+Train(:,i)*Train(:,i)'/L/tao_child(i);
        end
        R_MKA = abs(R_MKA_t);
        count =count+1;
        if count >=3  
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

