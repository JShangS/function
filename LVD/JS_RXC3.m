function [ tfr ] = JS_RXC3( x ,a )
%JS_RX3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%����x�ĶԳ�˲ʱ����غ��� ,a Ϊ�߶�ϵ��
L = length(x);
a = round(a);
%��ת
xf = x.';
xf = fliplr(xf);
% xf = xf.';
t =round( -L/2:L/2-1);
tao =round( -L/2:L/2-1);
% t = 0:L-1;
% LL = length(t);
L_tao = length(tao);
tfr = zeros(L_tao,L);
for itao = 1 : L_tao %�̶�taoȻ����tΪ�Ա���x(t+(tao+a)/2)x^H(t-(tao+a)/2)
    tao(itao);
    tao_t = round((tao(itao)+a/2));
    t1 = JS_BU0(x,tao_t,0);
    t2 = JS_BU0(x,-tao_t,0);
    tfr(itao,:) = (conj(t2).*(t1));
end
end

