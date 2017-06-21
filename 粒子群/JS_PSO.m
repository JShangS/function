clc
clear 
close all
%%%基本粒子群算法，根据《Particle Swarm Optimization》 James Kennedy
%%By Jerry Shang,JS
%%2017/6/5
N = 10;                                    %粒子群规模，即粒子数
D = 2;                                     %粒子的维数，即寻优参数的个数
X_pso = 1*rand(D,N);                     %初始化粒子群的位置
V = rand(D,N);                             %初始化粒子群的速度
Pbest = 10000*ones(1,N);                   %粒子群的个体最优解
Gbest = 0;                                 %粒子群的全体最优解索引
Max_iter = 100 ;                           %最大迭代次数
Xmax = 3;                                 %X的最大取值范围[-1,1]
Vmax = 10;                                 %最大速度
iter = 1;
g_increament = 0.01;                         %步长

tic
while(iter <= Max_iter)
    bool_Xmax = abs(X_pso) < Xmax;
    X_pso = X_pso.*bool_Xmax + (ones(D,N) - bool_Xmax) .* (Xmax * ones(D,N));
    Pbest_t = f2(X_pso(1,:),X_pso(2,:));
    bool_P = Pbest_t<Pbest;                %把Pbest中的大于当前X_pso结果的值替换掉
    Pbest =  Pbest_t.*bool_P + (ones(1,N) - bool_P) .* Pbest;
    Gvalue_iter(iter)=min(Pbest);           %每一代的最小值
    Gbest = find(Pbest==min(Pbest));        %找到当前的全局最优所在的位置
    X_iter(:,iter) = X_pso(:,Gbest);    %每一代的最好值
    bool_X = X_pso<X_pso(:,Gbest);  %向当前最好的X的位置靠拢
    direct = double(bool_X);
    direct = 2*(direct - 0.5);
    V =  V + direct.*rand(D,N)*g_increament;
    V(:,Gbest)=0;
%     for i = 1:D
%         bool_X = X_pso(i,:)<X_pso(i,Gbest);  %向当前最好的X的位置靠拢
%         direct = double(bool_X);
%         direct = 2*(direct - 0.5);
%         V(i,:) = V(i,:) + direct.*rand(1,N)*g_increament;
%         V(i,Gbest)=0;
%         bool_V(i,:) = V(i,:)<Vmax;
%         V(i,:) = V(i,:).*bool_V(i,:) + (ones(1,N) - bool_V(i,:)) .* (Vmax * ones(1,N));     
%     end
    X_pso = X_pso + V;
    iter = iter + 1;
end
time1 = toc;
Gbest = find(Pbest==min(Pbest)); 

disp(['PSO的最佳值: ',num2str(Pbest(:,Gbest))])
disp(['收敛的坐标:','(',num2str(X_iter(1,end)),',',num2str(X_iter(2,end)),')']);
disp(['PSO的耗时: ',num2str(time1)])


%每一代最优值
figure()
plot(Gvalue_iter);
xlabel('迭代次数')
ylabel('最小值')

% figure()
% plot(X_iter(1,:),'k');
% hold on
% plot(X_iter(2,:),'g');
% plot(1:Max_iter,ones(1,Max_iter),'r')
% legend('X1值','X2值','最价值')
% xlabel('迭代次数')
% ylabel('坐标值')
figure()
plot(X_iter(1,:),'k');
hold on
plot(X_iter(2,:),'g');
legend('X1值','X2值')
xlabel('迭代次数')
ylabel('坐标值')


step = 0.005 ;
L = -Xmax:step:Xmax;
[X1,X2]=meshgrid(L,L);
tic
Z = f2(X1,X2);
minZ = min(min(Z));
[hang,lie] = find(Z==minZ);
time2 = toc;
if step >=0.005
   figure()
%    mesh(X1,X2,Z);
   contour(X1,X2,Z)
   title('函数图像')
   hold on 
   plot(X_iter(1,:),X_iter(2,:),'b-*');
   plot(X_iter(1,end),X_iter(2,end),'r*','LineWidth',2,'MarkerSize',12);
   plot(L(hang),L(lie),'ko','LineWidth',2,'MarkerSize',12)
   legend('函数值等高线','PSO搜索路径','PSO最终值位置','穷举法的最好值位置')
%    plot3(X_iter(1,:),X_iter(2,:),Gvalue_iter,'r-*');
end
disp(['步长为',num2str(step),'下穷举的最佳值: ',num2str(minZ)])
disp(['最佳值的坐标:','(',num2str(L(hang)),' , ',num2str(L(lie)),')']);
disp(['穷举情况下的耗时: ',num2str(time2)])
