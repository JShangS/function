function [ X ] = JS_SFT2( x , a, h,Ts)
%JS_SFT �˴���ʾ�йش˺�����ժҪ
%%����
%%��Lv��s Distribution: Principle, Implementation,Properties, and Performance��
%%��ʽ(22)��д
%%%  2017.4.25
%%%By Shang Zheran����JerryShang
% x:�����ź�m*n��
% ��߶ȵ�DFT,��ÿ������߶�DFT
% X ÿһ���ǿ�ʱ��Ƶ��ÿ��һ����tao��
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

