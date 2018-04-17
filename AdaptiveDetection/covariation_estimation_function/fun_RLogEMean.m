function [ R ] = fun_RLogEMean( X )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%%%Log-E的均值协方差
[N,L] = size(X);
logm_R = zeros(N,N);
for i = 1:L
%     Ri = X(:,i) * X(:,i)';
    Ri = fun_Positive(X(:,i),3);
    logm_R = logm_R + fun_Logm(Ri);
end
logm_R = logm_R/L;
R = fun_Expm(logm_R);
end

