function [ Pd_mc,Th ] = fun_Detection( f,n,SNRout,PFA,N,rouR )
%FUN_DETECTION 此处显示有关此函数的摘要
%   此处显示详细说明
str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
irouR = inv(rouR);
SNRnum=10.^(SNRout/10);
MonteCarloPd=1e4;
MonteCarloPfa=100/PFA;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'/sqrt(N); %%%%%% 系统导向矢量
%%%%%正式开始%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%门限计算%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = waitbar(0,'Please wait...');
[fM,~] = size(f);
R_estimation = zeros(N,N,fM);
T = zeros(fM,MonteCarloPfa);
h = waitbar(0,'Please wait...');
for i = 1:MonteCarloPfa
    waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
%%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
    Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
   for i_R = 1: fM
       %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
       R_estimation(:,:,i_R) = f{i_R,1}(Train);
       %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
alpha=sqrt(SNRnum/abs(s'*irouR*s)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
h = waitbar(0,'Please wait...');
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
    for i=1:MonteCarloPd 
       %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%检测信号
        x0=alpha(m)*s+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
        for i_R = 1: fM
           %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
           R_estimation_t = f{i_R,1}(Train);
           %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

