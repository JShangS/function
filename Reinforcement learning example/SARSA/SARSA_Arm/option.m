%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                        %%
%%%                  Written by Matthew Kontz                              %%
%%%                  Walla Walla College                                   %%
%%%                  Edward F. Cross School of Engineering                 %%
%%%                  February 2001                                         %%
%%%                  Simulation of a planar three link robot.              %%																								%%%
%%%                                                                        %%
%%%      This function reconfigures the screen to the four different       %%
%%%      options that you can choose using the push buttons at the         %%
%%%      bottom of the display window.                                     %%
%%%                                                                        %%
%%%      This function is called by demobot.  To use, first execute        %%
%%%      demobot.  There are five files need to run demobot: demobot.m,    %%
%%%      option.m, forkin, invkin and setplot.                             %%
%%%                                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = option(OP,action)
global PushBut1 PushBut2 PushBut3 PushBut4
global tx txA Chose r psi rmax T1 T2 T3 rmin
global slider1 slider2 slider3 slider4 slider5
global r_last psi_last Down
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if OP<5 & action~=1
   set(PushBut1,'BackgroundColor',[0.8 0.8 0.8])
	set(PushBut2,'BackgroundColor',[0.8 0.8 0.8])	
	set(PushBut3,'BackgroundColor',[0.8 0.8 0.8])
	set(PushBut4,'BackgroundColor',[0.8 0.8 0.8])
	set(txA,'visible','off')
	set(slider1,'visible','off')
	set(slider2,'visible','off')
	set(slider3,'visible','off')  
	set(tx,'visible','off')
	set(slider4,'visible','off')
	set(slider5,'visible','off')
	set(gcf,'WindowButtonMotionFcn','')
	set(gcf,'WindowButtonDownFcn','')
   set(gcf,'WindowButtonUpFcn','')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OP==1    		% option 1 is position sliders
	r=rmax;
	psi=pi/2;
	invkin(r,psi)
	forkin
	setplot
	r=rmax;
	psi=pi/2;
	Chose=1;
	set(slider4,'value',rmax)
	set(slider5,'Value',pi/2)
	set(PushBut1,'BackgroundColor',[0.6 0.6 0.6])
	set(tx,'visible','on')
	set(slider4,'visible','on')
	set(slider5,'visible','on')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif OP==2		% option 2 is to click on target
  	% do option #2
	r=rmax;
	psi=pi/2;
	if action==0 
   	Chose=2;
		set(PushBut2,'BackgroundColor',[0.6 0.6 0.6])
		set(txA,'visible','off')
		% reset mouse button to nothing
   	set(gcf,'WindowButtonDownFcn','option(2,1)')
   	set(gcf,'WindowButtonMotionFcn','')
		set(gcf,'WindowButtonUpFcn','')
		invkin(r,psi)
		forkin
		setplot
   	set(gcf,'WindowButtonDownFcn','option(2,1)')
   	r_last=r;
		psi_last=psi;   
	else   
		step=10;			% number of steps between Tlast and T goal   
   	cp=get(gca,'CurrentPoint');
      set(gcf,'WindowButtonDownFcn','')
   	x=cp(1,1);
		y=cp(1,2);
		r=sqrt(x^2+y^2);
		psi=atan2(y,x);
   	if r<=rmax & r>=rmin & psi>=0 & psi<=pi
   		FLIP=sign((psi-pi/2)*(psi_last-pi/2));	%negative if robot flips
  			if FLIP==1;
  				Dr=(r-r_last)/step;
  				Dpsi=(psi-psi_last)/step;
  				for k=1:step;
  					invkin(r_last+k*Dr,psi_last+k*Dpsi)
  					forkin;
   	  			setplot
    	  			pause(0);
     			end  
 			else
  	  			T1_last=T1;
  				T2_last=T2;
  	  			T3_last=T3;
      		invkin(r,psi)
      		T1_goal=T1;
      		T2_goal=T2;
 				T3_goal=T3;
  				DeltaT1=(T1_goal-T1_last)/step;
  				DeltaT2=(T2_goal-T2_last)/step;
   			DeltaT3=(T3_goal-T3_last)/step;
  				for n=1:step;		% animated robot from last to goal
     				T1=T1_last+n*DeltaT1;
     				T2=T2_last+n*DeltaT2;
     				T3=T3_last+n*DeltaT3;
       	 		forkin
      	  		setplot
     				pause(0)
  				end
   		end      
   		r_last=r;
			psi_last=psi;
   	end
   	set(gcf,'WindowButtonDownFcn','option(2,1)')
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
elseif OP==3  			% this is for the click and drag option
	% reset plot
   T1=0;T2=0;T3=0;
   r_last=rmax;
   psi_last=pi/2;
	forkin
	setplot
	% reset figure window
	Chose=3;
	set(PushBut3,'BackgroundColor',[0.6 0.6 0.6])
	% do option #3
	set(gcf,'WindowButtonDownFcn','option(5,1);Down=1;')	
	set(gcf,'WindowButtonUpFcn','option(5,0);Down=0;')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif OP==4  			% this options is for angle sliders
	% resets plot
	r=rmax;
	psi=pi/2;
	invkin(r,psi)
	forkin
	setplot
	% reset figure window
	set(PushBut4,'BackgroundColor',[0.6 0.6 0.6])
	set(slider1,'Value',0)
	set(slider2,'Value',0)
	set(slider3,'Value',0)
	set(txA,'visible','on')
	set(slider1,'visible','on')
	set(slider2,'visible','on')
	set(slider3,'visible','on')
	Chose=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% this is the realtime inverse kinematics part of option #3
elseif OP==5				
	if action==1			%	Turns on click and drag option
	   set(gcf,'WindowButtonMotionFcn','option(5,666)')
	elseif action==0  	%  Turns off click and drag option
	   set(gcf,'WindowButtonMotionFcn','')
	else                 %  Makes robot tip follow mouse
      set(gcf,'WindowButtonMotionFcn','')					% turns of motion fcn
      set(gcf,'WindowButtonDownFcn','Down=1;')		% keep track of if Button up or down
		set(gcf,'WindowButtonUpFcn','Down=0;')		% and doesn't can any new functions	
      cp=get(gca,'CurrentPoint');
		x=cp(1,1);
		y=cp(1,2);
      r=sqrt(x^2+y^2);
		psi=atan2(y,x);
		if r<=rmax & r>=rmin & psi>=0 & psi<=pi   
         if (psi_last>=pi/2 & psi<pi/2) | (psi_last<pi/2 & psi>=pi/2)
            invkin(r_last,psi_last)
            step=6;
            T1_last=T1;
  				T2_last=T2;
 	  			T3_last=T3;
      		invkin(r,psi)
      		T1_goal=T1;
      		T2_goal=T2;
 				T3_goal=T3;
  				DeltaT1=(T1_goal-T1_last)/step;
  				DeltaT2=(T2_goal-T2_last)/step;
   			DeltaT3=(T3_goal-T3_last)/step;
  				for n=1:step;		% animated robot from last to goal
     				T1=T1_last+n*DeltaT1;
     				T2=T2_last+n*DeltaT2;
     				T3=T3_last+n*DeltaT3;
       	 		forkin
      	  		setplot
     				pause(0)
  				end
         else
            invkin(r,psi)
				forkin
				setplot
         	pause(0)	
         end
         r_last=r;
			psi_last=psi;
      end
      if Down==1
         set(gcf,'WindowButtonMotionFcn','option(5,666)') 
      else
         set(gcf,'WindowButtonMotionFcn','')
      end
      set(gcf,'WindowButtonDownFcn','option(5,1);Down=1;')	
		set(gcf,'WindowButtonUpFcn','option(5,0);Down=0;')
 	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   