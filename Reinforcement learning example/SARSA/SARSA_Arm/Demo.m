function Demo()
clc
clf;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global x1 y1 x2 y2 x3 y3 xt yt Pt			% position variable
global S1 S2 S3 S4 S5 S6					% position strings 
global pF1 pF2 pF3							% handles for fill
global dis Down          					% handles for text display
global C2 C3 Ct txA tx						% handles for the joint's circle
global J2 J3 Jt 							% handles for the joint's pluses
global L1 L2 L3 Link1 Link2 Link3			% Link matrices
global T1 T2 T3 STOP Chose					% input variables
global l1 l2 l3 rmax rmin Bmax				% constants
global TxtEpisode TxtSteps goal f1 f2 grafica
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = subplot(2,1,1);
f2 = subplot(2,1,2);
grafica = false;
subplot(f1);
P2 = ['setgrafica();'];
PushBut2=uicontrol(gcf,'Style','togglebutton','Units','normalized', ...
   	'Position',[0.8 .8 0.17125 0.05],'string','Step Graph', ...
      'Callback',P2,'visible','on','BackgroundColor',[0.8 0.8 0.8]);
set(gcf,'name','Reinforcement Learning with a three link planar robot');
set(gcf,'Color','w')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Base	
angle=225:15:315;		% get rit of where links overlap
xa=1.125*cos(angle*pi/180);
ya=1.125*sin(angle*pi/180);
xb=[xa .875 .875 2.5  2.5 -2.5 -2.5 -.875 -.875];
yb=[ya -.707107 -1.25 -1.25 -2.875 -2.875 -1.25 -1.25 -.707107];
B =[xb' yb' zeros(size(xb))' ones(size(xb))']';  % Link3 matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Link1 and Link2  
angle=265:-15:105;		% define an arc from 270 to 90 degree r=1.125
xa=1.125*cos(angle*pi/180);
ya=1.125*sin(angle*pi/180);
angle=150:15:210;		% get rit of where links overlap
xc=8.625+1.125*cos(angle*pi/180);
yc=1.125*sin(angle*pi/180);
xL1=[0 7.1875 7.1875 7.72978 xc 7.72978 7.1875 7.1875 0 xa 0];
yL1=[1.125 1.125 .875 .68133 yc -.68133 -.875 -1.125 -1.125 ya 1.125];
Link1=[xL1' yL1' zeros(size(xL1))' ones(size(xL1))']';  % Link1 matrix
Link2=Link1;				% Link1 and Link2 are the same
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Link3  
xL3=[0 2.875 3.125 3.125 4.875 6.625 6.6258 4.875 3.125 3.125 2.875 0 xa 0];
yL3=[1.125*ones(1,3) .875 .875 .25 -.25 -.875 -.875 -1.125*ones(1,3) ya 1.125];
Link3=[xL3' yL3' zeros(size(xL3))' ones(size(xL3))']';  % Link3 matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l1=8.625;				% distance between frame '1' and '2'
l2=l1;					% distance between frame '2' and '3'
l3=6.125;				% distance between frame '3' and 'tool'
rmax=l1+l2+l3;						% maximum distance between (0,0) and tool frame
rmin=(l2^2+(l1-l3)^2)^0.5;		% minimum distance between (0,0) and tool frame
Bmax=atan2(l1-l3,l2)+pi/2;
x0=0;						% x position of frame '0'	
y0=0;						% y position of frame '0'
x1=x0;					% x position of frame '1'	
y1=y0;					% y position of frame '1'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_i=rmax;				% initial radius
psi_i=90*pi/180;		% initial angle
T1=0;T2=0;T3=0;
forkin
r=rmax;
psi=pi/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xlabel('X-Position(in)')
%ylabel('Y-Position(in)')
axis([-25 25 -10 25])	    % axis limits
axis manual					% set axis to exact manual value(i.e [-25 25 -10 25])
axis equal					% x-scale=y-scale
hold on						% does not erase previous graphs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arc=(0:1:180)*pi/180;		% plot desired workspace
arc2=(180:-1:0)*pi/180;
plot([rmax*cos(arc) rmin*cos(arc2) rmax],[rmax*sin(arc) rmin*sin(arc2) 0], 'Color',[.8 .8 .8])
%legend('Workspace')
%Goal Position
goal = plot(0,0,'ok','markersize',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid off						% turns on grid
pos=[15,20];
%lg=legend('Workspace',1); % plot workspace or grid
%set(lg,'Position',[0.66 0.815 0.203571 0.0492857])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'XTick',[-25:5:25])	 % numbers on y-axis
set(gca,'YTick',[0:5:25])	 % numbers onf x-axis
%set(gca,'Color',[1,1,1]) 	 % plot background color
set(gca,'FontSize',7);				
%set(gca,'Color',[.95,.95,.95])	% edge background color
set(gca,'Position',[0.05 0.37 0.8 0.6])	% size of data windown
subplot(f2)
set(gca,'FontSize',7);
set(gca,'Position',[0.07 0.03 0.9 0.30])	% size of data windown
subplot(f1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BF=fill(B(1,:),B(2,:),'k');					        % color fill base
%set(BF,'FaceColor',[.8 .3 .3]);
set(BF,'FaceColor',[1 1 1]);
pF1=fill(L1(1,:),L1(2,:),'k', 'erasemode','normal');	% color fill link1
set(pF1,'FaceColor',[1 1 1]);
%set(pF1,'FaceColor',[.9 .9 .9]);
pF2=fill(L2(1,:),L2(2,:),'k', 'erasemode','normal');	% color fill link2
%set(pF2,'FaceColor',[.8 .9 .2]);
set(pF2,'FaceColor',[1 1 1]);
pF3=fill(L3(1,:),L3(2,:),'k','erasemode','normal');	% color fill link3
%set(pF3,'FaceColor',[.2 .2 .8]);
set(pF3,'FaceColor',[1 1 1]);
plot(x1,y1,'ok');					            % circle at joint '1'
C2=plot(x2,y2,'ok', 'erasemode','normal');		% circle at joint '2'
C3=plot(x3,y3,'ok', 'erasemode','normal');		% circle at joint '3'
Ct=plot(xt,yt,'ok', 'erasemode','normal');		% circle at joint 'T'
plot(x1,y1,'+k');					            % plus at joint '1'
J2=plot(x2,y2,'+k', 'erasemode','normal');		% plus at joint '2'
J3=plot(x3,y3,'+k', 'erasemode','normal');		% plus at joint '3'
Jt=plot(xt,yt,'+k', 'erasemode','normal');		% plus at joint 'T'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dis(1)=fill([-15.5 -15.5 15.5 15.5 -15.5],[-9 -4 -4 -9 -9],'w');    % plot a white box
dis(2)=plot([-15.5 -15.5 15.5 15.5 -15.5],[-9 -4 -4 -9 -9],'k');	% box's black outline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yrange = axis;
vspace = (yrange(4) - yrange(3))/20;
px=-14;
py=-5.5;
dis(3)=text(px, py,            	    mat2str(S1),'erasemode','normal');		% theta 1
dis(4)=text(px,(py-1.2*vspace), 	mat2str(S2),'erasemode','normal');		% xt
dis(5)=text(px+11, py,            	mat2str(S3),'erasemode','normal');		% theta2
dis(6)=text(px+11,(py-1.2*vspace),  mat2str(S4),'erasemode','normal');		% yt
dis(7)=text(px+21, py,            	mat2str(S5),'erasemode','normal');		% theta3
dis(8)=text(px+21,(py-1.05*vspace), mat2str(S6),'erasemode','normal');		% phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
px=-24.5;
py=-7;
TxtEpisode = text(px, py,'Episode: ','erasemode','normal','FontSize',8);			
TxtSteps   = text(px,(py-1*vspace),'Steps: ','erasemode','normal','FontSize',8);		

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gco,'BackingStore','off')					% for realtime inverse kinematics
set(gco,'Units','data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         REINFORCEMENT LEARNING LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow;

ArmDemo(200);


