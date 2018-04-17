function [ R ] = fun_RLogMeanPer( X )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%%%Log-E的均值协方差,基于反对称先验
[N,L] = size(X);
J = fun_permutation(N);
logmR = 0;
for i = 1:L
%     Ri = X(:,i) * X(:,i)';
    Ri = fun_Positive(X(:,i),3);
    logmR = logmR + fun_Logm(Ri + J * conj(Ri) * J);
end
% logmR = logm(X*X' + J * conj(X*X') * J) + logm(0.5*eye(N));

R = fun_Expm(logmR/L+logm(0.5*eye(N)));
end


