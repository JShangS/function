function [ R ] = fun_RLogEMean( X, opt)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%%%Log-E的均值协方差
[N,L] = size(X);
logm_R = zeros(N,N);
if nargin<2
    opt = 1;%%归一化的
end
if opt == 1
    for i = 1:L
    %     Ri = X(:,i) * X(:,i)';
        Ri = fun_Positive(X(:,i),4);
    %     logm_R = logm_R + fun_Logm(Ri);
        logm_R = logm_R + logm(Ri/(X(:,i)'*X(:,i)/N));
    end
else
    for i = 1:L
    %     Ri = X(:,i) * X(:,i)';
        Ri = fun_Positive(X(:,i),4);
    %     logm_R = logm_R + fun_Logm(Ri);
        logm_R = logm_R + logm(Ri);
    end
end

logm_R = logm_R/L;
% R = fun_Expm(logm_R);
R = expm(logm_R);
end

