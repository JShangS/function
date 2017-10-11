clc
clear 
close all
% �������ס�A New CFAR Detection Test for Radar����д��
% ����Ӧƥ���˲���AMF��
%%��ʱ��Ϊÿ������ĵ������ͬ��
%%%������MN*��ѹ���棬M=1��N=1,�źŷ���Ϊ1ʱ
%%%SNR����ȶ�Ӧ�ļ�����=SNR-10*log10��MN*��ѹ���棩ʱ�ļ�����
fc = 1e9;       %��Ƶ
C = 3e8;        %����
lamda = C/fc;   %����
tao = 127e-6;    %����
B = 1e6;        %����
mu = B/tao;     %��Ƶ��
Fs = 2*B;       %����Ƶ��
Ts = 1/Fs;
t = -tao/2:Ts:tao/2-Ts;  %��ʱ��
L = length(t);
R = 0e3;
dt = 2*R/C;
M = 2;                  %������
%%���в���
N = 1;                  %��Ԫ��
d = 0.5*lamda;          %��Ԫ���
theta = 0;%jiao2hu(10);    %�����
St = exp(1j*2*pi*d*(0:N-1).'/lamda*sin(theta));%����ʸ��
S = St;
S = repmat(S,[M,1]);    %����ʸ��
%%%���մ���(��ѹ)
signal = exp(-1j*2*pi*(fc*t+0.5*mu*t.^2) );
h_fft = fft(signal);
MC = 1000;
%%%
Pfa = 1e-6;     %�龯��
r0 = -log(Pfa); %����
SNR = -17:0.1:0;
for SNR_i = 1:length(SNR)
    SNR_i
    num = 0;
        for MC_i = 1:MC
            Pc = zeros(N,M*L);
%             A = 1*rand()+1j*1*rand();
            A=1;
            for i =1:M %����
                for j =1:N %����
                    echo = A*exp(1j*2*pi*(fc*dt+0.5*mu*(t-dt).^2) );
                    echo = awgn(echo,SNR(SNR_i),'measured');
                    echo_fft = fft(echo);
                    Pct = ifftshift(ifft(h_fft.*echo_fft)).*repmat(St(j),[1,L]);
                    Pct = Pct.*repmat(St(j),[1,L]);
                    Pc(j,1+(i-1)*L:i*L) = Pct;
                end
            end
          detect_R=L/2+1;
    %     detect_R=round(L/4);%%ûĿ���λ��
            x1 = [];
            x2 = [];
            x = [];
            y = [];
            for i =1:M
                x1 = Pc(:,1+L*(i-1):detect_R+L*(i-1)-10);
                x2=  Pc(:,detect_R+L*(i-1)+10:L*i);
                x = cat(1,x,cat(2,x1,x2));
                y = cat(1,y,Pc(:,detect_R+L*(i-1)));
            end
            M_gu = x*x'/length(x);
            r = abs((abs(S.'*inv(M_gu)*y))^2/(S.'*inv(M_gu)*S));%GLRT��AMF�����
            if r > r0
                num = num+1;
            end
        end
        proportion(SNR_i) = num/MC;
end
plot(SNR,proportion)

% figure()
% plot(abs(Pc(1,:)))
% %%%����ֵ
% %%%
Pfa = 1e-6;     %�龯��
r0 = -log(Pfa); %����
SNR = 0:0.1:20;
alpha = SNR2real(SNR,10);
n = 0:150;
for i =1:length(alpha)
    Pd(i) = sum(alpha(i).^n.*exp(-alpha(i))./factorial(n).*gammainc(n+1,r0));
end
% figure()
hold on
plot(SNR,Pd,'r')


