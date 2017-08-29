function [ LVD ] = JS_LVD(s,Ts,a,h)
%JS_LVD 此处显示有关此函数的摘要
%   此处显示详细说明
%%lv分布实现函数，经过瞬时自相光，变尺度FFT，在进行FFT后得到lv分布结果
%s:信号
%a:延迟，默认为1s
%h：默认为1
%Ts:采样频率
if (nargin<2)
    error('至少要输入信号和采样频率');
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

