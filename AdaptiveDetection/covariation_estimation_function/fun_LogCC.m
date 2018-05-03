function [ R_CC,alpha0 ] = fun_LogCC(X,R,R_KA )
%FUN_CC 此处显示有关此函数的摘要
%   此处显示详细说明
%%迭代求解，训练样本估计的协方差和先验协方差的线性组合，利用凸优化得到组合系数。
%%X:训练样本
%R,样本估计的协方差
%R_KA:先验协方差
% R = fun_SCMN(X);
R_CC = R;
[N,K] = size(X);
Rim = zeros(N,N,K);

for i = 1:K
    Rim(:,:,i) = logm(fun_Positive( X(:,i),4));
    
end
alpha = 0;
alpha0 =0;
for k = 1:20
    alpha = alpha0;
    t1 = 0;
    t2 = 0;
%     t3 = zeros(N,N);
    for i = 1:K 
        t1 = t1 + norm( Rim(:,:,i) - logm(R_CC),'fro')^2;
        t2 = t2 + norm(logm(R_KA) - Rim(:,:,i),'fro')^2;
    end
    alpha0 = max(min(1,(t1/t2)),0);
%     for i = 1:K
%         t3 = t3 + alpha0 * R_KA +(1-alpha0) * X(:,i) * X(:,i)';
%     end
%     R_CC_0 = R_CC;
%     R_CC = t3/K;    
    R_CC_0 = R_CC;
    R_CC = expm(alpha0 * logm(R_KA) +(1-alpha0) * logm(R_CC_0));
    if abs(alpha0- alpha)<1e-2
        break;
    end
end
end


