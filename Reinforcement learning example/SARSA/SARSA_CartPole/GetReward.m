function [ r,f ] = GetReward( s )
% r: the returned reward.
% f: true if the car reached the goal, otherwise f is false
    
x         = s(1);
x_dot     = s(3);
theta     = s(3);
theta_dot = s(4);

r = 10 - 10*abs(10*theta)^2 - 5*abs(x) - 10*theta_dot;
f = false;


twelve_degrees     = deg2rad(12); % 12º
fourtyfive_degrees = deg2rad(45); % 45º
%if (x < -4.0 | x > 4.0  | theta < -twelve_degrees | theta > twelve_degrees)          
if (x < -4.0 | x > 4.0  | theta < -fourtyfive_degrees | theta > fourtyfive_degrees)          
    r = -10000 - 50*abs(x) - 100*abs(theta);     
    f = true;
end

    
   


    
