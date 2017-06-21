function AcrobotPlot( x,a,steps )
subplot(2,1,2);

theta1 = x(1);
theta2 = x(2);


x_acrobot(1)=0;
y_acrobot(1)=0;

x_acrobot(2) = x_acrobot(1) + sin(theta1); 
y_acrobot(2) = y_acrobot(1) - cos(theta1);

x_acrobot(3) = x_acrobot(2) + sin(theta2); 
y_acrobot(3) = y_acrobot(2) - cos(theta2); 
	
plot(x_acrobot ,y_acrobot,'ok-','LineWidth',1,'markersize',7,'MarkerFaceColor',[.7 .7 .7]);
hold on
plot(x_acrobot(3),y_acrobot(3),'.r','markersize',20);

title(strcat ('Step: ',int2str(steps)));

%-----------------------
axis([-2.1 2.1 -2.1 2.1])
axis square
grid on
drawnow
hold off