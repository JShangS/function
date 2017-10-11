%opt2d.m: ȫ���ɶȿ�ʱ����Ӧ����
%--------------------------------------------------------------------------
%start         : 2004.11.04  AM 10:21   yunrisheng
%Latest change : 
%--------------------------------------------------------------------------
clc
close all
clear
tic
%�Ӳ��������
N = 12;                        % ��Ԫ����
M = 10;                        % ���������
CNR = 30;                      % �����
beta = 1;                      % �Ӳ��۵�ϵ��(beta = 2*v*T/d)
sita_a = -90:.9:90.;           % �Ӳ���Ԫ����               
sita = sita_a*pi/180;
[NN N_bin] = size(sita);
%Ŀ�����
sita_t = -25;                  % Ŀ��DOA
omiga_t = 0.4;                 % Ŀ��Doppler
SNR = -20;                       % �����

%�ռ䵼��ʸ����ʱ�䵼��ʸ��
%�ռ�Ƶ�ʺ�DoppleƵ������ omiga_d = beta * omiga_s
omiga_s = pi*sin(sita);      
omiga_d = beta*omiga_s;       

aN = zeros(N,N_bin);
bN = zeros(M,N_bin);

aN = exp(-j*[0:N-1]'*omiga_s)./sqrt(N); %%��λ
bN = exp(-j*[0:M-1]'*omiga_d)./sqrt(M); %%Doppler


%Ŀ���ʱ�ź�
aN_t = zeros(N,1);
bN_t = zeros(M,1);

aN_t = exp(-j*pi*[0:N-1]'*sin(sita_t*pi/180))/sqrt(N);
bN_t = exp(-j*pi*[0:M-1]'*omiga_t)/sqrt(M);

S_t = zeros(M*N,1);
S_t = kron(aN_t,bN_t);

%�����Ӳ�Э�������
R = zeros(M*N,M*N);                     
S = zeros(M*N,N_bin);                   
ksai = 10^(CNR/10)*(randn(1,N_bin)+j*randn(1,N_bin))/sqrt(2);               %������̬�ֲ��������ֵ������Ϊ1
for ii = 1:N_bin
    S(:,ii) = kron(aN(:,ii),bN(:,ii));  
    R = R + ksai(ii).*(S(:,ii)*S(:,ii)');       
end

%����Э������������Ϊ30dB
R = R +eye(M*N);     %CNR = 30dB
inv_R = inv(R);                   %�����
%������ֵ��
[u s v] = svd(R);                       %����ֵ�ֽ�
figure(1);
plot(10*log10(diag(s)));
title('��Ԫ��N=12, ���������M=10');
axis([0 120 -10 50]);
xlabel('����ֵ��Ŀ');
ylabel('����ֵ(dB)');
grid on

P_f = zeros(N_bin,N_bin);
P_min_var = zeros(N_bin,N_bin);
%���Ӳ���
for ii = 1:N_bin
    for jj = 1:N_bin
            SS = kron(aN(:,ii),bN(:,jj));
            P_f(ii,jj) = SS'*R*SS;        %������
            P_min_var(ii,jj) = 1./(SS'*inv_R*SS);
    end
end        

%��С�������
figure(2)
mesh(sin(sita),omiga_d/pi,20*log10(abs(P_min_var)));
title('��Ԫ��N=12, ���������M=10');
xlabel('��λ����');
ylabel('��һ��DoppleƵ��');
zlabel('����(dB)');
grid on
%��ʱ����Ȩ����
tic
w_opt = inv(R)*S_t./(S_t'*inv_R*S_t);
%w_opt = inv(RR)*a_t;

%�����ſ�ʱ��Ӧ
for ii = 1:N_bin
    for jj = 1:N_bin
        SSS = kron(aN(:,ii),bN(:,jj));
        res_opt(ii,jj) = SSS'*w_opt;
    end
end
t_stap=toc
                      
figure(3)
%[X,Y]=meshgrid(omiga_d/pi,sita_a);
% imagesc(omiga_d/pi,omiga_d/pi,abs(res_opt))
% mesh(omiga_d/pi,omiga_d/pi,abs(res_opt))
mesh(omiga_d/pi,omiga_d/pi,10*log10(abs(res_opt).^2))
title('��Ԫ��N=12, ���������M=10');
xlabel('��һ��oppleƵ��');
ylabel('��λ����');
zlabel('����(dB)');
zlabel('����');
grid on
%�����Ÿ�������
for ii = 1:N_bin
    for jj = 1:N_bin
        SS = kron(aN(:,ii),bN(:,jj));
 %       IF(ii,jj) = SS'*inv_R*SS.*trace(R)./(SS'*SS);
        IF(ii,jj) = SS'*inv_R*SS./(SS'*SS);
    end
end
figure(4)
% axis([-1 1 10 55]);
% mesh(10*log10(abs(IF)));
plot(omiga_d/pi,10*log10(abs(IF(101,:))));
xlabel('��һ��oppleƵ��');
ylabel('��������(dB)');
grid on



%%%AMF
tic
for ii = 1:N_bin
    for jj = 1:N_bin
        SSS = kron(aN(:,ii),bN(:,jj));
        res_opt_AMF(ii,jj) = abs(abs(S_t'*inv_R*SSS)^2/(S_t'*inv_R*S_t));
    end
end
t_AMF=toc
figure(5)
% mesh(omiga_d/pi,omiga_d/pi,abs(res_opt_AMF))
% imagesc(omiga_d/pi,omiga_d/pi,abs(res_opt_AMF))
% contour(omiga_d/pi,omiga_d/pi,abs(res_opt_AMF),5000)
mesh(omiga_d/pi,omiga_d/pi,10*log10(abs(res_opt_AMF).^2))
title('��Ԫ��N=12, ���������M=10');
xlabel('��һ��oppleƵ��');
ylabel('��λ����');
zlabel('����(dB)');
% zlabel('����');
grid on
[hang,lie] = find(max(max(abs(res_opt_AMF)))==abs(res_opt_AMF));
figure(6)
plot(omiga_d/pi,10*log10(abs(res_opt_AMF(:,lie)).^2))
