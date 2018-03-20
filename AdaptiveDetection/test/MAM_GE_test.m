%%%多先验模型和信息几何度量实现测试2018.3.20
clc
clear
close all
rou0=0.945;
% rou=0.5:0.05:0.99;
rou=[0.945,0.9,0.1];
N=8;
L=round(2*N);
LL=1000;
R0=fun_rho(rou0,N);
for i =1:length(rou)
    R=fun_rho(rou(i),N);
    MAM(:,:,i)=R;
end
x0 = fun_TrainData('p',N,1,R0,3,1,1);
Rx0 = x0*x0';
[V,D] = svd(Rx0);
Evalue = abs(diag(D));
index_1 = find(Evalue<=1);
Evalue(index_1) = 1;
% Evalue = sort(Evalue,'descend');
D = diag(Evalue);
Rx0 = (abs(V)*D*abs(V'));
[R_MAM_r,~,ratio_r]=fun_information_estimation(Rx0,MAM,'c')
