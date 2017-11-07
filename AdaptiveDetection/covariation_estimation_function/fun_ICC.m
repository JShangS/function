function [ R_ICC ] = fun_ICC( X,R2,R1 )
%FUN_ICC 此处显示有关此函数的摘要
%   此处显示详细说明
%%迭代的CC
R1 = fun_CC(X,R2,R1);
end

