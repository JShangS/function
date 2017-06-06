function [ CRB_result ] = JS_CRB_PolyPhase( I, N, Ts, SNR )
%CRB_ 此处显示有关此函数的摘要
% % 多项式相位信号的CRB
% % I：阶次,0:初相、1:中心频率、2:调频率、....
% % N：估计时的信号点数
% % Ts： 采样间隔
% % SNR：信噪比
% % CRB_result： 每一行是不同信噪比下的单个参数估计下限
D = diag (Ts.^(0:I));
n0 = -N/2 + 1/2;
n =n0 : N/2 - 1/2;
H = zeros(I,I);
for k = 0:I
    for l = 0:I
        H(k+1,l+1) = sum(n.^(k+l));
    end
end
SNR = 10.^(-SNR./10);
inv_D = inv(D);
inv_H = inv(H);
CRB_result = zeros(I+1, length(SNR));
for i = 1:length(SNR)
    CRB_result(:,i) = diag(0.5 * SNR(i) * inv_D * inv_H * inv_D);
end
end

