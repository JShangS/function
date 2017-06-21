%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                        %%
%%%                  Written by Matthew Kontz                              %%
%%%                  Walla Walla College                                   %%
%%%                  Edward F. Cross School of Engineering                 %%
%%%                  February 2001                                         %%
%%%                  Simulation of a planar three link robot.              %%																								%%%
%%%                                                                        %%
%%%      The purpose of this function is to calculate all the values       %%		
%%%      needed to update the plot for a given theta1, theta2 and theta3.  %%
%%%                                                                        %%
%%%      This function is called by demobot.  To use, first execute        %%
%%%      demobot.  There are five files need to run demobot: demobot.m,    %%
%%%      option.m, forkin, invkin and setplot.                             %%
%%%                                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = forkin()
global x1 y1 x2 y2 x3 y3 xt yt Pt			% position variable
global S1 S2 S3 S4 S5 S6					% position strings 
global L1 L2 L3 Link1 Link2 Link3			% Link matrices
global T1 T2 T3 									% input variables
global l1 l2 l3 									% constants6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P1=pi/2+T1*pi/180;	% angle between matlab and frame '1'
P2=P1+T2*pi/180;		% angle between matlab and frame '2'
P3=P2+T3*pi/180;		% angle between matlab and frame '3'
Pt=P3;					% angle between matlab and frame 'tool'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x2=l1*cos(P1);			% x position of frame '2'	
y2=l1*sin(P1);			% y position of frame '2'
x3=x2+l2*cos(P2);		% x position of frame '3'	
y3=y2+l2*sin(P2);		% y position of frame '3'
xt=x3+l3*cos(P3);		% x position of frame 'tool'	
yt=y3+l3*sin(P3);		% y position of frame 'tool'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tm1 = [ 	cos(P1) 	-sin(P1)		0		x1		% Transform from matlab to frame '1'
   		sin(P1)  cos(P1) 		0		y1
   			0			0			1		0
         	0			0			0		1 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tm2 = [ 	cos(P2) 	-sin(P2)		0		x2 	% Transform from matlab to frame '2'
   		sin(P2)  cos(P2) 		0		y2
   			0			0			1		0
         	0			0			0		1 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tm3 = [ 	cos(P3) 	-sin(P3)		0		x3		% Transform from matlab to frame '3'
   		sin(P3)  cos(P3) 		0		y3
   			0			0			1		0
         	0			0			0		1 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
L1=Tm1*Link1;
L2=Tm2*Link2;
L3=Tm3*Link3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S1=['\theta_{1} = ' num2str(0.01*round(100*T1),3) '^o'];
S2=['X_{T} =' num2str(0.01*round(100*xt),3) ''];
S3=['\theta_{2} = ' num2str(0.01*round(100*T2),3) '^o'];
S4=['Y_{T} = ' num2str(0.01*round(100*yt),3) ''];
S5=['\theta_{3} = ' num2str(0.01*round(100*T3),3) '^o'];
S6=['\phi_{T} = ' num2str(0.01*round(100*Pt*180/pi),3) '^o'];
