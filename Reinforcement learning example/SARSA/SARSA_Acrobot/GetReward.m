function [ r,f ] = GetReward( x )
% MountainCarGetReward returns the reward at the current state
% x: a vector of position and velocity of the car
% r: the returned reward.
% f: true if the car reached the goal, otherwise f is false


theta1 = x(1);
theta2 = x(2);
y_acrobot(1) = 0;
y_acrobot(2) = y_acrobot(1) - cos(theta1);
y_acrobot(3) = y_acrobot(2) - cos(theta2);    




%goal
goal = y_acrobot(1) + 1.0 ;

r = -1;%y_acrobot(3);
f = false;

if( y_acrobot(3) >= goal) 
	r = 100;%10*y_acrobot(3);        
    f = true;
end

   


    
   


    
