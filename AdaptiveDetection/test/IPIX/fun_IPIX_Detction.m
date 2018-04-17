function [ Pd_mc,Th ] = fun_IPIX_Detction( f,n,SNRout,PFA)
%FUN_IPIX_DETCTION 此处显示有关此函数的摘要
%   此处显示详细说明
%%%f:用到的协方差估计器，检测器，和总体名称
%%%是个M（检测器个数）行，N（协方差，检测器，姓名）列的元孢数组
% clc
% clear 
% close all
Read_Display_Data
Data_process
load(matFile) 
%%%%参数设置
SNRnum=10.^(SNRout/10);
MonteCarloPd=1e4;
MonteCarloPfa=100/PFA;
rouR = R_KA;  %%真实的杂波协方差
irouR = inv(rouR);
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'/sqrt(N); %%%%%% 系统导向矢量
%%%%%正式开始%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%门限计算%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = waitbar(0,'Please wait...');
Zhh = sig;
tic
[fM,~] = size(f);
R_estimation = zeros(N,N,fM);
T = zeros(fM,MonteCarloPfa);
h = waitbar(0,'Please wait...');
for i = 1:MonteCarloPfa
    waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
%%%%%%%%用实测数据确定门限，但是数据有限。%%%%%%%%%%%%%%%%%%%%%%%
    index_t1 = ceil(rand()*(M-10));
    Train1 = Zhh(index_t1:index_t1+7,Range-L/2:Range-1);
    Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2);
    Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = Zhh(index_t1:index_t1+N-1,Range) ; % 接收信号仅包括杂波和噪声
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        index_t1 = ceil(rand()*(M-10));
        Train1 = Zhh(index_t1:index_t1+7,Range-L/2:Range-1);
        Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2);
        Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        x0 = Zhh(index_t1:index_t1+7,Range) ; % 接收信号仅包括杂波和噪声
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
