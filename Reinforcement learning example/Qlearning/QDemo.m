clc
clear 
close all

Rew=zeros(20);Rew(19,15)=100;    % rewards
task=struct('initialState',[1,1],'terminalState',[19,15]); % task
robot=struct('alpha',1,'gamma',1,'Qtable',zeros(400,4),'best',[],'state',[1,1]);
robot=Qlearning(robot,task,Rew,200);

disp('前驱状态 | 动作 | 后继状态')
for l=1:size(robot.best,1)
    s0(l,:)=robot.best(l,[1,2]);
    a(l,:)=robot.best(l,3);
    s(l,:)=robot.best(l,[4,5]);
    disp(['  ',num2str(s0(l,:)),'   |  ',num2str(a(l,:)),'   |  ',num2str(s(l,:))]);
end
contour(Rew)
hold on
plot(s(:,2),s(:,1),'r')