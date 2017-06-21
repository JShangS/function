]clear 
close all
%%%��������Ⱥ�㷨�ĸĽ� Current Simplified particle swarn optimizer
%%%�����ݡ�Particle Swarm Optimization�� James Kennedy
%%% ����Ҫ��g_increament �����������ã��������û��һ���Ľ�����ʽ�����ԱȽϿ����顣
%%It was soon realized that there is no good way to guess whether p - or g-increment should be larger. 
%%By Jerry Shang,JS
%%2017/6/8
N = 100;                                %����Ⱥ��ģ����������
D = 2;                                 %���ӵ�ά������Ѱ�Ų����ĸ���
prsentx = 1-2*rand(D,N);                 %��ʼ������Ⱥ��λ��,������Ϊ��һ����λ��
pbestx = rand(D,N);                      %���ֵ��Ӧ��λ�� 
V = rand(D,N);                             %��ʼ������Ⱥ���ٶ�
Pbest = 10000*ones(1,N);                    %����Ⱥ�ĸ������Ž�
Gbest = 0;                                 %����Ⱥ��ȫ�����Ž�����
Max_iter = 100 ;                           %����������
Record_X = zeros(D,N,Max_iter);              %��¼ÿһ�������ӷֲ�
Xmax = 2;                                  %X�����ȡֵ��Χ[-1,1]
Vmax = 10;                                 %����ٶ�
w = 0.9;                                   %����Ȩֵ
c1 = 1;                                    %�������� 
c2 = 1;
iter = 1;
X_iter = zeros(N,Max_iter);
Gvalue_iter = zeros(1,Max_iter);
tic
while(iter <= Max_iter)
    Record_X(:,:,iter) = prsentx;
    %%XԽ����ֹ
    bool_Xmax = abs(prsentx) < Xmax;
    prsentx = prsentx.*bool_Xmax + (ones(D,N) - bool_Xmax) .* (Xmax * ones(D,N));
    %%%
    Pbest_t = f1(prsentx(1,:),prsentx(2,:));  %��һ����ֵ
    bool_P = Pbest_t<Pbest;                 %��Pbest�еĴ��ڵ�ǰpbestx�����ֵ�滻��
    %����ÿ��������õ�ֵ
    Pbest =  Pbest_t.*bool_P + (ones(1,N) - bool_P) .* Pbest;
    %����ÿ��������õ�λ��
    pbestx =  prsentx.*repmat(bool_P,2,1) + repmat((ones(1,N) - bool_P),2,1) .* pbestx;
    Gvalue_iter(iter)=min(Pbest);           %ÿһ������Сֵ
    Gbest = find(Pbest==min(Pbest));        %�ҵ���ǰ��ȫ���������ڵ�λ��
    Gt =  pbestx(:,Gbest);
    [Gt_hang,Gt_lie] = size(Gt);
    Gt = reshape(Gt,Gt_hang*Gt_lie,1);
    X_iter(1:Gt_hang*Gt_lie,iter) = Gt;       %ÿһ�������ֵ
    V = w*V + c1*rand(D,N).*(pbestx - prsentx) + c2*rand(D,N).*(pbestx(:,Gbest) - prsentx);
    V(:,Gbest)=0;
    %%����V�Ĵ�С
%     bool_V = abs(V)<Vmax;
%     V = V.*bool_V + (ones(D,N) - bool_V) .* (Vmax * ones(D,N));
    %��һ����x
    prsentx = prsentx + V; 
    %������1
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

figure()
plot(X_iter(1,:),'k');
hold on
plot(X_iter(2,:),'g');
legend('X1ֵ','X2ֵ')
xlabel('��������')
ylabel('����ֵ')


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


% figure()
% for i = 1:Max_iter
%     plot(Record_X(1,:,i),Record_X(2,:,i),'r.')
%     hold on
%     axis([-10,10,-10,10])
%     pause(0.001)
%     
% end
