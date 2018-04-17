function [ R_CC,alpha0 ] = fun_CCIter(X,R_KA )
%FUN_CC 此处显示有关此函数的摘要
%   此处显示详细说明
%%迭代求解，训练样本估计的协方差和先验协方差的线性组合，利用凸优化得到组合系数。
%%X:训练样本
%R,样本估计的协方差
%R_KA:先验协方差
R = fun_SCMN(X);
[R_CC,alpha0_1] = fun_CC(X,R_KA);
while(1)
    [R_CC,alpha0] = fun_CC2(X,R_CC,R_KA);
    if abs(alpha0_1 - alpha0) <1e-5
        break;
    end
    alpha0_1 = alpha0;
end
% iter = 20;
% alpha0 = zeros(iter+1,1);
% [R_CC,alpha0(1)] = fun_CC(X,R,R_KA);
% for i = 1 :iter
%     [R_CC,alpha0(i+1)] = fun_CC(X,R_CC,R_KA);
% end
end


