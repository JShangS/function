function [ y ] = Q( x )
%   高斯右尾函数
%   x为积分下限，y为积分值，标准正态分布X>x的概率 
y = 0.5*erfc(x/sqrt(2));
end

