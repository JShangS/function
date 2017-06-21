function  ArmDemo( maxepisodes )
%ArmDemo, the main function of the demo
%maxepisodes: maximum number of episodes to run the demo
global TxtEpisode goal f1 f2 grafica MOVIE T1 T2 T3

%MOVIE = getframe(gcf);

clc
statelist   = BuildStateList();  % the list of states
actionlist  = BuildActionList(); % the list of actions

nstates     = size(statelist,1);
nactions    = size(actionlist,1);

Q1           = BuildQTable( nstates,nactions );  % the QTable
Q2           = BuildQTable( nstates,nactions );  % the QTable
Q3           = BuildQTable( nstates,nactions );  % the QTable

% no learning (stable system)
%x =load('Qtables.mat');
%Q1 = x.Q1;
%Q2 = x.Q2;
%Q3 = x.Q3;
%alpha = 0.0;
%epsilon     = 0.01;

maxsteps    = 400;  % maximum number of steps per episode
alpha       = 0.3;  % learning rate
gamma       = 1.0;  % discount factor
epsilon     = 0.01; % probability of a random action selection
grafica     = false; % indicates if display the graphical interface


xs=0;
ys=23.4;
xpoints=[];
ypoints=[];
xf = 10;
yf = 15;
%[xf yf]= randgoal();
set(goal,'xdata',xf,'ydata',yf);

goals = randgoalArray();
N = size(goals,1);
goal_index=1;
cum_r=0;

T1 = 0;
T2 = 0;
T3 = 0;

for i=1:maxepisodes 
    
%     if (goal_index<=N)
%         xf = goals(goal_index,1);
%         yf = goals(goal_index,2);    
%         goal_index = goal_index+1;
%         goal_index=1; %cyclic trainning
%     else
         [xf yf]= randgoal();  % activate the random location of the goal        
%     end
    
    absedist = edist([xf yf],[xs ys]);
    
    
    set(goal,'xdata',xf,'ydata',yf);
    
    
    
    set(TxtEpisode,'string',strcat('Episode: ',int2str(i)));     
    [total_reward,steps,Q1,Q2,Q3 ] = ArmEpisode( maxsteps, Q1, Q2, Q3, xf , yf, alpha, gamma,epsilon,statelist,actionlist,grafica );    
    
    if (mod(i,20)==0)
        save Qtables.mat Q1 Q2 Q3;
    end
    
    disp(['Espisode: ',int2str(i),' steps: ',int2str(steps),' reward: ',num2str(total_reward),' eplison:',num2str(epsilon)])
    
    xpoints(i) = i-1;
    %ypoints(i) = steps/absedist; %attempt to ignore the distance from the
    %starting point.
    ypoints(i) = steps;
    %cum_r = cum_r + total_reward;
    %ypoints(i)  = cum_r;
    subplot(f2);
    plot(xpoints,ypoints);
    xlabel('Episodes');
    ylabel('total reward');
    %legend('Learning curve','Location','NorthWest');
    subplot(f1);
    %text(xf,yf,int2str(i));
    %plot(xf,yf,'k');
    setplot;
    drawnow;
    epsilon = epsilon * 0.99;
    
    if (i>1000000)
        grafica=true;
    end
    %MOVIE(numel(MOVIE)+1) = getframe(gcf);
    
end

%movie2avi(MOVIE,'Arm.avi');




