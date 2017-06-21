function [ s ] = DiscretizeState( x, statelist )
%DiscretizeState check which entry in the state list is more close to x and
%return the index of that entry.

x(1) = sign(x(1));
x(2) = sign(x(2));
x(4) = sign(x(4));

x = repmat(x,size(statelist,1),1);
[d  s] = min(edist(statelist,x));

%[d  s] = min(dist(statelist,x'));
