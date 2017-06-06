function [ X ] = JS_SFT( x , Ts , fai)
%JS_SFT �˴���ʾ�йش˺�����ժҪ
%%����
%%��Lv��s Distribution: Principle, Implementation,Properties, and Performance��
%%��ʽ(22)��д
%%%  2017.4.25

% x:�����ź�m*n��
% Ts����ʵ�Ĳ������ 
% a�� h���߶�ϵ����Ĭ��1
% ��߶ȵ�DFT,��ÿһ������߶�DFT
% X ÿһ���ǿ�ʱ��Ƶ��ÿ��һ�б�ĳ߶�
%%
if(nargin < 2)
    error('������������������źźͲ���ʱ����');
end
if(nargin == 2)
    fai = 1;
end
if(nargin == 3)
    h = 1;
end
[N, M] = size ( x ); %�źų���N,M���ӳ�
X = zeros(N, M);
for im = 1 : M
    for in = 1 : N
        exp_t = exp(-1j * 2 * pi * fai * in / N * Ts * (0:N-1)*Ts);
        X(in , im) = exp_t * x(:,in);
    end
end
X = ifft(X,[],1);
end

