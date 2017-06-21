function [ r,f ] = ArmGetReward( x, xf, yf )
% GetReward returns the reward at the current state
% x: a vector of position of the last arm of the robot
% r: the returned reward.
% f: true if the car reached the goal, otherwise f is false
  

penalty = 1.0;
f       = false;
T1      = x(1);
T2      = x(2);
T3      = x(3);
xt      = x(4);
yt      = x(5);

dist = (xf-xt)^2 + (yf-yt)^2;


%  if ( xf>0 && T2>0  )
%       penalty =  penalty + abs(T2);
%  end
%   
%  if ( xf>0 && T3>0 )
%       penalty =  penalty + abs(T3);
%  end
% 
% if ( xf<0 && T2<0   )
%      penalty =  penalty + abs(T2);
% end
% 
% if ( xf<0 && T3<0 )
%      penalty =  penalty + abs(T3);
% end

  if ( xf>0 && sign(T2)>0 && sign(T3)>0  )
       penalty =  10 + abs(T2) + abs(T3);
  end
  
  if ( xf<0 && sign(T2)<0 && sign(T3)<0  )
       penalty =  10 + abs(T2) + abs(T3);
  end


r = -((dist^3)/100.0) + (penalty);
r = r /10;
%r = -1 - penalty;


if ( sqrt(dist)<1 && penalty==1)    
    r =  (10^13)/(1+dist);         
    f =  true;
    %r = 1000;
end

    
   


    
