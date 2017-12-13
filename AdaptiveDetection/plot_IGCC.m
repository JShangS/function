%%%Äægamma·Ö²¼µÄÍ¼
clc
clear 
close all
N = 1;
L = 1000;
M = 1;
lambda = linspace(1,5,5);
mu = 1;
opt_train = 1;
figure(1)
hold on
for i = 1:length(lambda)
    h = ( abs(fun_TrainData_IGCC( N,L,M,lambda(i),mu,opt_train)));
    [f,xi] = ksdensity(h);
    [n,~] = hist(h,100);
    n= n/L;
    if i == 1
        plot(xi,f,'k','linewidth',2);
    elseif i== 2
        plot(xi,f,'g','linewidth',2);
    elseif i== 3
        plot(xi,f,'b','linewidth',2);
    elseif i== 4
        plot(xi,f,'c','linewidth',2);
    elseif i== 5
        plot(xi,f,'y','linewidth',2);
    end
end
% axis([1,100,0,1])
h_leg = legend('\lambda=1','\lambda=2','\lambda=3','\lambda=4','\lambda=5');
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
% set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
