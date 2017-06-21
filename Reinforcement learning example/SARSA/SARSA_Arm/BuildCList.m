function [ CList ] = BuildCList(states,actions )
%BuildCList Build the classifiers list
%states: the list of states features
%actions: the initial action list


Ns = size(states,1);
Nf
Na = size(actions,1);

index=1;
CList =[];
for i=1:Ns
    for j=1:Na
        CList(index,
        

