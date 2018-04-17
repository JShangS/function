function [ Pd_mc,Th ] = fun_Detection( f,n,SNRout,PFA,N,rouR )
%FUN_DETECTION �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
str_train = 'p';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
irouR = inv(rouR);
SNRnum=10.^(SNRout/10);
MonteCarloPd=1e4;
MonteCarloPfa=100/PFA;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'/sqrt(N); %%%%%% ϵͳ����ʸ��
%%%%%��ʽ��ʼ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%���޼���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = waitbar(0,'Please wait...');
[fM,~] = size(f);
R_estimation = zeros(N,N,fM);
T = zeros(fM,MonteCarloPfa);
h = waitbar(0,'Please wait...');
for i = 1:MonteCarloPfa
    waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
%%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
    Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % �����źŽ������Ӳ�������
   for i_R = 1: fM
       %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%
       R_estimation(:,:,i_R) = f{i_R,1}(Train);
       %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       T(i_R,i) = f{i_R,2}(R_estimation(:,:,i_R),x0,s);
   end   
end
close (h)
T_descend=sort(T,'descend');
for i = 1 : fM
  Th(i) = (T_descend(i,floor(MonteCarloPfa*PFA-1))+T_descend(i,floor(MonteCarloPfa*PFA)))/2;  
end
counter = zeros(1,fM);
Pd_mc = zeros(fM,length(SNRout));
alpha=sqrt(SNRnum/abs(s'*irouR*s)); % ����SNR=|alpha|^2*s'*R^(-1)*s���|alpha|
h = waitbar(0,'Please wait...');
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
    for i=1:MonteCarloPd 
       %%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
        x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % �����źŽ������Ӳ�������
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%����ź�
        x0=alpha(m)*s+x0;%+pp;    %%%%%%%  ��Ҫ  %%%%%%%%%%%%%
        for i_R = 1: fM
           %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%
           R_estimation_t = f{i_R,1}(Train);
           %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           Tt = f{i_R,2}(R_estimation_t,x0,s);
           if Tt > Th(i_R)
              counter(i_R) =  counter(i_R) + 1;
           end
        end 
    end
    for i_Pd = 1 : fM
        Pd_mc(i_Pd,m) = counter(i_Pd)/MonteCarloPd; counter(i_Pd) = 0;
    end
end
close (h)
end

