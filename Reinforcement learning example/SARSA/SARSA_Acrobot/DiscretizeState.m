function [ s ] = DiscretizeState( x, statelist )
%DiscretizeState check which entry in the state list is more close to x and
%return the index of that entry.
%x(3) = sign(x(3));
%x(4) = sign(x(4));

%[d  s] = min(dist(statelist,x'));

x = repmat(x,size(statelist,1),1);
[d  s] = min(edist(statelist,x));