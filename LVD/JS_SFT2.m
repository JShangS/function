function [ X ] = JS_SFT2( x , a, h,Ts)
%JS_SFT 此处显示有关此函数的摘要
%%根据
%%《Lv’s Distribution: Principle, Implementation,Properties, and Performance》
%%公式(22)编写
%%%  2017.4.25
%%%By Shang Zheran——JerryShang
% x:输入信号m*n，
% 变尺度的DFT,对每行做变尺度DFT
% X 每一行是快时间频域，每是一列随tao变
%%%[(tao1,t1),(tao1,t2),(tao1,t3)
%%% (tao2,t1),(tao2,t2),(tao2,t3)
%%% (tao3,t1),(tao3,t2),(tao3,t3)]
%%
[M,N] = size ( x );
m = -M/2:M/2-1; %tao
n = -N/2:N/2-1; %t
l = m;%0:N-1;
X = zeros(M,N);
% for im = 1 : M
%     for in = 1 : N
%         exp_t = exp(-1j * 2 * pi * (l(im)) * n * (2*m(in)*Ts + a)*h/ N);%.* exp(-1j * 2 * pi * 2 *l(im) * n * (m(in)*Ts + a) * h/ N) ;
%         X(im , in) = x(in,:) * exp_t.';
%     end
% end
t1 = n.'*(2*m*Ts + a)*h/N;
for im =1:M
    im
    exp_t = exp(-1j*2*pi*l(im)*t1);
    X(im,:) = diag(x * exp_t);
end


% X = (ifft(X,[],2));%ifftshift
end

