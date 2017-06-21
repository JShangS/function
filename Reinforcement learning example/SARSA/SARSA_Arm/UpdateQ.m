function [ Q ] = UpdateQ( s1, s2, r, a, tab , alpha, gamma )
% UpdateQ update de Qtable and return it using Whatkins QLearing
% s1: previous state before taking action (a)
% s2: current state after action (a)
% r: reward received from the environment after taking action (a) in state
%                                             s1 and reaching the state s2
% a:  the last executed action
% tab: the current Qtable
% alpha: learning rate
% gamma: discount factor
% Q: the resulting Qtable

Q = tab;
Q(s1,a) =  Q(s1,a) + alpha * ( r + gamma* max(Q(s2,:)) - Q(s1,a) );