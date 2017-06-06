function [ tfr ] = JS_RXC3( x ,a )
%JS_RX3 此处显示有关此函数的摘要
%   此处显示详细说明
%%计算x的对称瞬时自相关函数 ,a 为尺度系数
L = length(x);
a = round(a);
%翻转
xf = x.';
xf = fliplr(xf);
% xf = xf.';
t =round( -L/2:L/2-1);
tao =round( -L/2:L/2-1);
% t = 0:L-1;
% LL = length(t);
L_tao = length(tao);
tfr = zeros(L_tao,L);
for itao = 1 : L_tao %固定tao然后以t为自变量x(t+(tao+a)/2)x^H(t-(tao+a)/2)
    tao(itao);
    tao_t = round((tao(itao)+a/2));
    t1 = JS_BU0(x,tao_t,0);
    t2 = JS_BU0(x,-tao_t,0);
    tfr(itao,:) = (conj(t2).*(t1));
end
end

