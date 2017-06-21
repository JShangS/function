function plot_Cart_Pole(s,action,steps)

x     = s(1);
theta = s(3);
l=3;
%l=3/(1+abs(theta));     %pole's Length for ploting it can be different from the actual length

pxg = [x+1 x-1 x-1 x+1 x+1];
pyg = [0.25 0.25 1.25 1.25 0.25];

pxp=[x x+l*sin(theta)];
pyp=[1.25 1.25+l*cos(theta)];

arrowfactor_x=sign(action)*2.5;
if (sign(arrowfactor_x)>0)
    text_arrow = strcat('==>> ',int2str(10*action));
else if (sign(arrowfactor_x)<0)
        text_arrow = strcat(int2str(10*action),' <<==');
    else
        text_arrow='=0=';
        arrowfactor_x=0.25;
    end
end

    
    
subplot(2,1,1);

%Car 
fill(pxg,pyg,[.6 .6 .5],'LineWidth',2);  %car
hold on
title(['Steps: ',int2str(steps)])

%Car Wheels
plot(x-0.5,0.25,'rO','LineWidth',2,'Markersize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5]);
plot(x+0.5,0.25,'rO','LineWidth',2,'Markersize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5]);

%Pendulum
%plot(pxp,pyp,'-rO','LineWidth',5,'Markersize',10,'MarkerEdgeColor','r','MarkerFaceColor','r');
plot(pxp,pyp,'-r','LineWidth',5);
plot(pxp(1),pyp(1),'.r','LineWidth',2,'Markersize',10,'MarkerEdgeColor','k','MarkerFaceColor','r');
plot(pxp(2),pyp(2),'rO','LineWidth',2,'Markersize',15,'MarkerEdgeColor','k','MarkerFaceColor','r');

text(x + arrowfactor_x - 0.5 ,0.8,text_arrow);
%axis([x-6 x+6 0 6])
axis([-6 6 0 6])
%grid
box off
drawnow;
hold off
