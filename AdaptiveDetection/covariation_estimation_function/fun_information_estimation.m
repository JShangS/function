function [ R_MAM, distance,ratio] = fun_information_estimation(R0, MAM)
%FUN_INFORMATION_ESTIMATION 此处显示有关此函数的摘要
%   此处显示详细说明
%%%%%%基于信息几何度量的协方差估计方法。
%%%R0:参考协方差
%%%MAM:多模型组合协方差为N*N，L为模型个数
[N,~,L]=size(MAM); 
distance = zeros(L,1);
for i=1:L
    distance(i) = fun_ReimanDistance(R0,MAM(:,:,i));
end
distance = distance/sum(distance);
Sum_ratio = sum(exp(-distance));
ratio = exp(-distance)/Sum_ratio;
ratio = reshape(ratio,1,1,L);
ratio = repmat(ratio,N,N,1);
R_MAM = sum(ratio.*MAM,3);
% ratio = ratio(1,1,1:L);
end

