clc
clear 
close all
%%%��������Ⱥ�㷨�����ݡ�Particle Swarm Optimization�� James Kennedy
%%By Jerry Shang,JS
%%2017/6/5
N = 10;                                    %����Ⱥ��ģ����������
D = 2;                                     %���ӵ�ά������Ѱ�Ų����ĸ���
X_pso = 1*rand(D,N);                     %��ʼ������Ⱥ��λ��
V = rand(D,N);                             %��ʼ������Ⱥ���ٶ�
Pbest = 10000*ones(1,N);                   %����Ⱥ�ĸ������Ž�
Gbest = 0;                                 %����Ⱥ��ȫ�����Ž�����
Max_iter = 100 ;                           %����������
Xmax = 3;                                 %X�����ȡֵ��Χ[-1,1]
Vmax = 10;                                 %����ٶ�
iter = 1;
g_increament = 0.01;                         %����

tic
while(iter <= Max_iter)
    bool_Xmax = abs(X_pso) < Xmax;
    X_pso = X_pso.*bool_Xmax + (ones(D,N) - bool_Xmax) .* (Xmax * ones(D,N));
    Pbest_t = f2(X_pso(1,:),X_pso(2,:));
    bool_P = Pbest_t<Pbest;                %��Pbest�еĴ��ڵ�ǰX_pso�����ֵ�滻��
    Pbest =  Pbest_t.*bool_P + (ones(1,N) - bool_P) .* Pbest;
    Gvalue_iter(iter)=min(Pbest);           %ÿһ������Сֵ
    Gbest = find(Pbest==min(Pbest));        %�ҵ���ǰ��ȫ���������ڵ�λ��
    X_iter(:,iter) = X_pso(:,Gbest);    %ÿһ�������ֵ
    bool_X = X_pso<X_pso(:,Gbest);  %��ǰ��õ�X��λ�ÿ�£
    direct = double(bool_X);
    direct = 2*(direct - 0.5);
    V =  V + direct.*rand(D,N)*g_increament;
    V(:,Gbest)=0;
%     for i = 1:D
%         bool_X = X_pso(i,:)<X_pso(i,Gbest);  %��ǰ��õ�X��λ�ÿ�£
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

disp(['PSO�����ֵ: ',num2str(Pbest(:,Gbest))])
disp(['����������:','(',num2str(X_iter(1,end)),',',num2str(X_iter(2,end)),')']);
disp(['PSO�ĺ�ʱ: ',num2str(time1)])


%ÿһ������ֵ
figure()
plot(Gvalue_iter);
xlabel('��������')
ylabel('��Сֵ')

% figure()
% plot(X_iter(1,:),'k');
% hold on
% plot(X_iter(2,:),'g');
% plot(1:Max_iter,ones(1,Max_iter),'r')
% legend('X1ֵ','X2ֵ','���ֵ')
% xlabel('��������')
% ylabel('����ֵ')
figure()
plot(X_iter(1,:),'k');
hold on
plot(X_iter(2,:),'g');
legend('X1ֵ','X2ֵ')
xlabel('��������')
ylabel('����ֵ')


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
   title('����ͼ��')
   hold on 
   plot(X_iter(1,:),X_iter(2,:),'b-*');
   plot(X_iter(1,end),X_iter(2,end),'r*','LineWidth',2,'MarkerSize',12);
   plot(L(hang),L(lie),'ko','LineWidth',2,'MarkerSize',12)
   legend('����ֵ�ȸ���','PSO����·��','PSO����ֵλ��','��ٷ������ֵλ��')
%    plot3(X_iter(1,:),X_iter(2,:),Gvalue_iter,'r-*');
end
disp(['����Ϊ',num2str(step),'����ٵ����ֵ: ',num2str(minZ)])
disp(['���ֵ������:','(',num2str(L(hang)),' , ',num2str(L(lie)),')']);
disp(['�������µĺ�ʱ: ',num2str(time2)])
