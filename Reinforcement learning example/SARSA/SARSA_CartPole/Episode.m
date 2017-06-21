function [ total_reward,steps,Q ] = Episode( maxsteps, Q , alpha, gamma,epsilon,statelist,actionlist,grafic )


global grafica
% state variables x,x_dot,theta,theta_dot
x            = [0 0 0 0.01];
steps        = 0;
total_reward = 0;


% convert the continous state variables to an index of the statelist
s   = DiscretizeState(x,statelist);
%s    = get_box(x);
% selects an action using the epsilon greedy selection strategy
a   = e_greedy_selection(Q,s,epsilon);


for i=1:maxsteps    
        
    % convert the index of the action into an action value
    action = actionlist(a);    
    
    %do the selected action and get the next car state    
    xp  = DoAction( action , x );    
    
    % observe the reward at state xp and the final state flag
    [r,f]   = GetReward(xp);
    total_reward = total_reward + r;
    
    % convert the continous state variables in [xp] to an index of the statelist    
    sp  = DiscretizeState(xp,statelist);
    %sp  = get_box(xp);
    % select action prime
    ap = e_greedy_selection(Q,sp,epsilon);
    
    
    % Update the Qtable, that is,  learn from the experience
    Q = UpdateSARSA( s, a, r, sp, ap, Q , alpha, gamma );
    
    
    %update the current variables
    s = sp;
    a = ap;
    x = xp;
    
    
    %increment the step counter.
    steps=steps+1;
    
   grafic=grafica;
    % Plot of the mountain car problem
    if (grafic==true)        
       plot_Cart_Pole(x,action,steps);    
    end
    
    % if the car reachs the goal breaks the episode
    if (f==true)
        break
    end
    
end


