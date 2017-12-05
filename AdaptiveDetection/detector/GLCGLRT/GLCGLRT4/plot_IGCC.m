clc
clear 
close all
N = 1;
L = 10000000;
M = 1;
lambda = 1:0.2:2;
mu = 1;
opt_train = 1;
figure(1)
hold on
for i = 1:length(lambda)
    h = ( abs(fun_TrainData_IGCC( N,L,M,lambda(i),mu,opt_train)));
    [n,~] = hist(h,1000);
    n= n/L;
    if i == 1
        plot(n,'k','linewidth',2);
    elseif i== 2
        plot(n,'g','linewidth',2);
    elseif i== 3
        plot(n,'b','linewidth',2);
    elseif i== 4
        plot(n,'c','linewidth',2);
    elseif i== 5
        plot(n,'y','linewidth',2);
    end
end
axis([1,100,0,1])
h_leg = legend('\lambda=1','\lambda=2','\lambda=3','\lambda=4','\lambda=5');
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
