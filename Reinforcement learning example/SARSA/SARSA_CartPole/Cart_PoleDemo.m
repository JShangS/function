function  Cart_PoleDemo( maxepisodes )
%MountainCarDemo, the main function of the demo
%maxepisodes: maximum number of episodes to run the demo

% Cart Pole control demo with SARSA 
% Programmed in Matlab 
% by:
%  Jose Antonio Martin H. <jamartinh@fdi.ucm.es>

global grafica



maxsteps    = 1000;              % maximum number of steps per episode
statelist   = BuildStateList();  % the list of states
actionlist  = BuildActionList(); % the list of actions

nstates     = size(statelist,1);
nactions    = size(actionlist,1);
Q           = BuildQTable( nstates,nactions );  % the Qtable

alpha       = 0.3;   % learning rate
gamma       = 1.0;   % discount factor
epsilon     = 0.001;  % probability of a random action selection
grafica     = false; % indicates if display the graphical interface

xpoints     = [];
ypoints     = [];

for i=1:maxepisodes    
    
    [total_reward,steps,Q ] = Episode( maxsteps, Q , alpha, gamma,epsilon,statelist,actionlist,grafica ); 
    
    disp(['Espisode: ',int2str(i),'  Steps:',int2str(steps),'  Reward:',num2str(total_reward),' epsilon: ',num2str(epsilon)])
    
    epsilon = epsilon * 0.99;
    
    xpoints(i)=i-1;
    ypoints(i)=steps;    
    subplot(2,1,2);    
    plot(xpoints,ypoints)      
    title(['Episode: ',int2str(i),' epsilon: ',num2str(epsilon)])  
    xlabel('Episodes')
    ylabel('Steps')    
    drawnow
    
    if (i>1000)
        grafica=true;
    end
end






