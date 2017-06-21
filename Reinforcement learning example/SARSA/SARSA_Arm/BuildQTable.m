function [ Q ] = BuildQTable( nstates,nactions )
%BuildQTable do exactly this
%Q: the returned initialized QTable

%Q = zeros(nstates,nactions);
%Q = zeros(nstates,nactions)-1; % a variant
%Q = zeros(nstates,nactions)+100; % another variant
%Q =  100*rand(nstates,nactions); % a variant
Q = random('Normal',0,100,nstates,nactions); % a normal distribution variant