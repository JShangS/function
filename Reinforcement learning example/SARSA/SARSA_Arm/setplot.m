%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                        %%
%%%                  Written by Matthew Kontz                              %%
%%%                  Walla Walla College                                   %%
%%%                  Edward F. Cross School of Engineering                 %%
%%%                  February 2001                                         %%
%%%                  Simulation of a planar three link robot.              %%																								%%%
%%%                                                                        %%
%%%      The purpose of this function is to update the figure window       %%	
%%%      with the values, matrix and strings created by forkin.            %%
%%%                                                                        %%
%%%      This function is called by demobot.  To use, first execute        %%
%%%      demobot.  There are five files need to run demobot: demobot.m,    %%
%%%      option.m, forkin, invkin and setplot.                             %%
%%%                                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = setplot()
global x2 y2 x3 y3 xt yt 						% position variable
global S1 S2 S3 S4 S5 S6						% position strings 
global pF1 pF2 pF3								% handles for fill
global C2 C3 Ct									% handles for the joint's circle
global J2 J3 Jt									% handles for the joint's pluses
global dis											% handles for text display
global L1 L2 L3 Link1 Link2 Link3			% Link matrices
global STOP Chose T1 T2 T3
global f1 f2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(f1);
set(pF1,'xdata',L1(1,:),'ydata',L1(2,:));	% 3 fills (2-t)
set(pF2,'xdata',L2(1,:),'ydata',L2(2,:));    
set(pF3,'xdata',L3(1,:),'ydata',L3(2,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(C2,'xdata',x2,'ydata',y2);  		% 3 circles (2-t)
set(C3,'xdata',x3,'ydata',y3);   
set(Ct,'xdata',xt,'ydata',yt);  
set(J2,'xdata',x2,'ydata',y2);   		% 3 plus (2-t)
set(J3,'xdata',x3,'ydata',y3);   
set(Jt,'xdata',xt,'ydata',yt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(dis(3),'string',mat2str(S1)); 				% theta 1
set(dis(4),'string',mat2str(S2));				% xt
set(dis(5),'string',mat2str(S3));				% theta 2
set(dis(6),'string',mat2str(S4));				% yt
set(dis(7),'string',mat2str(S5));				% theta 3
set(dis(8),'string',mat2str(S6)); 				% phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axis([-25 25 -10 25])
end
