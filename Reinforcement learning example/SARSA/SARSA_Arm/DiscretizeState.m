function [ s ] = DiscretizeState( x, statelist  )
%DiscretizeState check which entry in the state list is more close to x and
%return the index of that entry.

x     = sign(x);
%x(1) = sign(x(1)); %difference in x
%x(2) = sign(x(2)); %difference in y
[d  s] = min(dist(statelist,x'));
