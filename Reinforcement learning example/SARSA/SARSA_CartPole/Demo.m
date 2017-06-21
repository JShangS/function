function Demo()
clc
clf;
clear all;

global TxtEpisode TxtSteps goal f1 f2 grafica
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = subplot(2,1,1);
box off

f2 = subplot(2,1,2);
grafica = false;


P2 = ['setgrafica();'];
PushBut2=uicontrol(gcf,'Style','togglebutton','Units','normalized', ...
   	'Position',[0.83 .9 0.16 0.05],'string','Graficar', ...
      'Callback',P2,'visible','on','BackgroundColor',[0.8 0.8 0.8]);
set(gcf,'name','Reinforcement Learning with a Scara Manipulator Robot');
set(gcf,'Color','w')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid off						% turns on grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gco,'BackingStore','off')  % for realtime inverse kinematics
set(gco,'Units','data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         REINFORCEMENT LEARNING LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow;

Cart_PoleDemo(300);