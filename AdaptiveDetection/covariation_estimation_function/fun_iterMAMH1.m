function [ X0 ] = fun_iterMAMH1( MAM,Z)
%FUN_ITERMAMH0 此处显示有关此函数的摘要
%   此处显示详细说明
%%%《Adaptive detection of distributed targets in compound-Gaussian clutter 
%%without secondary data: An approach based on multiple a-priori spectral
%%models》，H1条件下的估计
%%MAM:多模型
%%Z：检测单元数据
iter_num = 50;%%迭代次数
[N,L] = size(Z);
[~,~,num] = size(MAM);
X0 = zeros(N,N); %%%%初始协方差
for i = 1:num
    X0 = X0 + 1/num * MAM(:,:,i);%%% 初始协方差为多模型的均值
end
tao = zeros(L,1);
T = zeros(N,N,L);
t = zeros(num,1);
for iter = 1:iter_num
    X0_t = zeros(N,N);
    for i = 1:L
        tao(i) = Z(:,i)' * Z(:,i) / L;
        for j = 1:L
             T(:,:,i) = T(:,:,i) + Z(:,j) * Z(:,j)' / tao(i);
        end
        t(i) = N / trace(X0 * T(:,:,i) / L);
        X0_t = X0_t + t(i) * MAM(:,:,i);
    end 
    X0 = X0_t / L;
end
end

