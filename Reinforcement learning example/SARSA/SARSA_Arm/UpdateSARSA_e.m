function [Q, e] = UpdateSARSA_e( s, a, r, sp, ap, tab, etrace, alpha, gamma, lambda )
% UpdateQ update de Qtable and return it using SARSA
% s1: previous state before taking action (a)
% s2: current state after action (a)
% r: reward received from the environment after taking action (a) in state
%                                             s1 and reaching the state s2
% a:  the last executed action
% tab: the current Qtable
% alpha: learning rate
% gamma: discount factor
% Q: the resulting Qtable
% lambda: elegibility trace decay factor


Q = tab;
e = etrace;

delta  =  r + gamma * Q(sp,ap) - Q(s,a);
e(s,a) =  e(s,a) + 1;

Q = Q + alpha * delta * e;
e = gamma * lambda * e;


