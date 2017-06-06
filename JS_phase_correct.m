function [ x ] = JS_phase_correct( s )
% 此处显示有关此函数的摘要
%   此处显示详细说明
%%ISAR中相位矫正用的函数,非相参
%%s为矫正后的数据，每一列为是一次回波。
[N, M] = size(s);  %s的行数和列数,距离单元数和回波数
e = zeros(1,M);  %%矫正项
e(1) = 1;  %以第一列为基准
for i = 2 : M %每一列做差, 得到初相差 
    sum1 = sum(s(:, i) .* conj(s(:, i-1)));
    t = sum1 / abs(sum1);
    e(i) = e(i-1) * t;
end 
%%加窗
w = 0.54 - 0.46 * cos(2 * pi * (1:M) ./ M); %%对每一行加窗
ww = repmat(w, N, 1);  %扩充窗函数
ee = repmat(e, N, 1); %%扩充矫正项
x = s  .* ww .* conj(ee);
end

