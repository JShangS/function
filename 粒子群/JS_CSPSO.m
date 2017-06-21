]clear 
close all
%%%基本粒子群算法的改进 Current Simplified particle swarn optimizer
%%%，根据《Particle Swarm Optimization》 James Kennedy
%%% 不需要对g_increament 参数进行设置，这个参数没有一定的解析方式，所以比较靠经验。
%%It was soon realized that there is no good way to guess whether p - or g-increment should be larger. 
%%By Jerry Shang,JS
%%2017/6/8
N = 100;                                %粒子群规模，即粒子数
D = 2;                                 %粒子的维数，即寻优参数的个数
prsentx = 1-2*rand(D,N);                 %初始化粒子群的位置,迭代作为上一代的位置
pbestx = rand(D,N);                      %最好值对应的位置 
V = rand(D,N);                             %初始化粒子群的速度
Pbest = 10000*ones(1,N);                    %粒子群的个体最优解
Gbest = 0;                                 %粒子群的全体最优解索引
Max_iter = 100 ;                           %最大迭代次数
Record_X = zeros(D,N,Max_iter);              %记录每一代的粒子分布
Xmax = 2;                                  %X的最大取值范围[-1,1]
Vmax = 10;                                 %最大速度
w = 0.9;                                   %惯性权值
c1 = 1;                                    %加速因子 
c2 = 1;
iter = 1;
X_iter = zeros(N,Max_iter);
Gvalue_iter = zeros(1,Max_iter);
tic
while(iter <= Max_iter)
    Record_X(:,:,iter) = prsentx;
    %%X越界阻止
    bool_Xmax = abs(prsentx) < Xmax;
    prsentx = prsentx.*bool_Xmax + (ones(D,N) - bool_Xmax) .* (Xmax * ones(D,N));
    %%%
    Pbest_t = f1(prsentx(1,:),prsentx(2,:));  %下一代的值
    bool_P = Pbest_t<Pbest;                 %把Pbest中的大于当前pbestx结果的值替换掉
    %更新每个粒子最好的值
    Pbest =  Pbest_t.*bool_P + (ones(1,N) - bool_P) .* Pbest;
    %更新每个粒子最好的位置
    pbestx =  prsentx.*repmat(bool_P,2,1) + repmat((ones(1,N) - bool_P),2,1) .* pbestx;
    Gvalue_iter(iter)=min(Pbest);           %每一代的最小值
    Gbest = find(Pbest==min(Pbest));        %找到当前的全局最优所在的位置
    Gt =  pbestx(:,Gbest);
    [Gt_hang,Gt_lie] = size(Gt);
    Gt = reshape(Gt,Gt_hang*Gt_lie,1);
    X_iter(1:Gt_hang*Gt_lie,iter) = Gt;       %每一代的最好值
    V = w*V + c1*rand(D,N).*(pbestx - prsentx) + c2*rand(D,N).*(pbestx(:,Gbest) - prsentx);
    V(:,Gbest)=0;
    %%限制V的大小
%     bool_V = abs(V)<Vmax;
%     V = V.*bool_V + (ones(D,N) - bool_V) .* (Vmax * ones(D,N));
    %下一代的x
    prsentx = prsentx + V; 
    %迭代加1
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

figure()
plot(X_iter(1,:),'k');
hold on
plot(X_iter(2,:),'g');
legend('X1值','X2值')
xlabel('迭代次数')
ylabel('坐标值')


step = 0.01 ;
L = -Xmax:step:Xmax;
[X1,X2]=meshgrid(L,L);
tic
Z = f1(X1,X2);
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


% figure()
% for i = 1:Max_iter
%     plot(Record_X(1,:,i),Record_X(2,:,i),'r.')
%     hold on
%     axis([-10,10,-10,10])
%     pause(0.001)
%     
% end
