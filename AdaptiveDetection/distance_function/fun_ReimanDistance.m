function [ distance ] = fun_ReimanDistance( A,B )
%FUN_REIMANDISTANCE 此处显示有关此函数的摘要
%   此处显示详细说明
%%%两个矩阵的黎曼几何距离，根据文献<基于信息几何的雷达目标检测_于方浦>公式4.18和4.20编的
T = A^(-0.5)*B*A^(-0.5);
[V,D] = eig(T);%%%特征分解，V为特征向量，D为特征值对角阵。
lnD = abs(log(D).^2);
distance = trace(lnD);
end

