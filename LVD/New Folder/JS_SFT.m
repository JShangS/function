function [ X ] = JS_SFT( x , Ts , fai)
%JS_SFT 此处显示有关此函数的摘要
%%根据
%%《Lv’s Distribution: Principle, Implementation,Properties, and Performance》
%%公式(22)编写
%%%  2017.4.25

% x:输入信号m*n，
% Ts：真实的采样间隔 
% a， h：尺度系数，默认1
% 变尺度的DFT,对每一列做变尺度DFT
% X 每一列是快时间频域，每是一行变的尺度
%%
if(nargin < 2)
    error('输入参数至少两个，信号和采样时间间隔');
end
if(nargin == 2)
    fai = 1;
end
if(nargin == 3)
    h = 1;
end
[N, M] = size ( x ); %信号长度N,M个延迟
X = zeros(N, M);
for im = 1 : M
    for in = 1 : N
        exp_t = exp(-1j * 2 * pi * fai * in / N * Ts * (0:N-1)*Ts);
        X(in , im) = exp_t * x(:,in);
    end
end
X = ifft(X,[],1);
end

