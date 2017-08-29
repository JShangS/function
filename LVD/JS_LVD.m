function [ LVD ] = JS_LVD(s,Ts,a,h)
%JS_LVD �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%lv�ֲ�ʵ�ֺ���������˲ʱ����⣬��߶�FFT���ڽ���FFT��õ�lv�ֲ����
%s:�ź�
%a:�ӳ٣�Ĭ��Ϊ1s
%h��Ĭ��Ϊ1
%Ts:����Ƶ��
if (nargin<2)
    error('����Ҫ�����źźͲ���Ƶ��');
end
if (nargin==2)
    a = 1;
    h = 1;
end
q = a/Ts;
RXC = JS_RXC3(s,q);
S = JS_SFT2(RXC,a,h,Ts); 
LVD = fftshift(fft(S,[],2),2);
end

