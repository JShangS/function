function [ s ] = DiscretizeState( x, statelist )
%DiscretizeState check which entry in the state list is more close to x and
%return the index of that entry.



%[d  s] = min(dist(statelist,x'));

x = repmat(x,size(statelist,1),1);
[d  s] = min(edist(statelist,x));

