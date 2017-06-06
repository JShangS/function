function [ tfr ] = JS_RXC2( x ,a)
%JS_RXC 此处显示有关此函数的摘要
%   此处显示详细说明
if (nargin == 1)
    a = 0;
end
M = length(x);
tfr = zeros(M, M);
tao = -M/2 : M/2;
for it = 1:M
    for itao = 1:length(tao)
        t_tao = round(tao(itao)/2 + a/2);
        if (it-t_tao)>0 & (it-t_tao)<=M &(it+t_tao)>0 & (it+t_tao)<=M 
            tfr(itao,it) = conj(x(it - t_tao)) * x(it + t_tao);
        end
    end
    t_in=round(1) ;
    if (it - t_in) >= 1 &(it + t_in) <= M & (it + t_in)>0 &(it - t_in) <=M
    tfr(t_in,it) = 0.5 * ( conj(x(it - t_in)) * x(it + t_in) + ...
        conj(x(it + t_in)) * x(it - t_in) );
    end
end
% tfr = fftshift(tfr,1);
end

