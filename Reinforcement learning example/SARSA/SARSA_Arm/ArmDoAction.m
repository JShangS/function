function [ xp ] = ArmDoAction( force )
%MountainCarDoAction: executes the action (a) into the mountain car
%environment
% a: is the force to be applied to the car
% x: is the vector containning the position and speed of the car
% xp: is the vector containing the new position and velocity of the car

global T1 T2 T3 xt yt

%bounds for angles
maxangle =  90;
minangle = -90;


T1 = T1 + force(1);
T2 = T2 + force(2);
T3 = T3 + force(3);

x =[ T1 T2 T3];

ind = find(x>=90);
x(ind)=90;

ind = find(x<=-90);
x(ind)=-90;

T1 = x(1);
T2 = x(2);
T3 = x(3);



forkin;

xp=[];
xp(1) = T1;
xp(2) = T2;
xp(3) = T3;
xp(4) = xt;
xp(5) = yt;






