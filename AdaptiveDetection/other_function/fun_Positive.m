function [ Rx ] = fun_Positive( X,opt )
%FUN_POSITIVE 此处显示有关此函数的摘要
%   此处显示详细说明
% 当矩阵不正定时化成正定的。
if nargin == 1
    opt = 1;
end
[N,~] = size(X);
if opt == 1
%%《Covariance matrix estimation via geometric
%%barycenters and its application to radar training data selection》
    Rx = X*X';
    [V,D] = svd(Rx);
    Evalue = abs(diag(D));
    index_1 = find(Evalue<=1);
    Evalue(index_1) = 1;
    % Evalue = sort(Evalue,'descend');
    D = diag(Evalue);
    Rx = abs(V)*D*abs(V');
elseif opt == 2
%%《Geometric barycenters for covariance estimation in compound-Gaussian
%%clutter》
    Rx = X*X';
    [V,D] = svd(Rx);
    Evalue = abs(diag(D));
    Km = max(Evalue)/min(Evalue);
    xk = norm(X,'fro');
    lambdak = max(1, xk^2*Km/(Km^2+N-1) );
    index_1 = find(Evalue<=lambdak);
    Evalue(index_1) = lambdak;
%     Evalue(1) = lambdak * Km;
    D = diag(Evalue);
    Rx = abs(V)*D*abs(V');
    Rx = Rx./xk^2;
elseif opt == 3
   Rx = (X*X'+ eye(N)); 
elseif opt ==4
   Rx = fun_corrcoef(X);
end
end

