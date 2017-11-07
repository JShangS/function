function [ R_AIWCM ] = fun_AIWCM( X,x0 )
%FUN_AIWCM 此处显示有关此函数的摘要
%   此处显示详细说明
% 《Adaptively iterative weighting covariance matrix estimation 
% for airborne radar clutter suppression》，2015
%X：训练样本
%x0：检测单元样本
[M,N] = size(X);
R_AIWCM = eye(M,M);%R_x0;%eye(N,N);
count = 1;
beta_child = ones(N,1);
beta_parent = zeros(N,1);
while (norm(beta_parent-beta_child,'fro')/norm(beta_child,'fro')>0.1)
    beta_parent = beta_child;
    R_AIWCM_inv = inv(R_AIWCM);
    beta_child_t1 = X'*R_AIWCM_inv*x0;
    beta_child_t2 = diag(X'*R_AIWCM_inv*X);
    beta_child = abs(beta_child_t1./beta_child_t2);
    R_AIWCM_t = 0;
    for i = 1:N
        R_AIWCM_t = R_AIWCM_t+beta_child(i)*(X(:,i)*X(:,i)')/sum(beta_child);%sum(beta_child);
    end
    R_AIWCM = (R_AIWCM_t);
    count =count+1;
    if count >50
        break;
    end
end
end

