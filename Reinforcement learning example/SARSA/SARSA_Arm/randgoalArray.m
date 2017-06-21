function [ goals ] = randgoalArray()

l1=8.625;				            % distance between frame '1' and '2'
l2=l1;					            % distance between frame '2' and '3'
l3=6.125;				            % distance between frame '3' and 'tool'
rmax=l1+l2+l3;						% maximum distance between (0,0) and tool frame
rmin=(l2^2+(l1-l3)^2)^0.5;		    % minimum distance between (0,0) and tool frame

goals=[];
index=1;
signo = 1;
for i=-1:0.05:1
    for j=-1*signo:0.05*signo:1*signo
        phi = (i)*pi;
        radius = ((rmax-rmin)/2.0)*(j)+((rmax+rmin)/2.0);
        goals(index,1)= abs(radius * cos(phi));
        goals(index,2)= abs(radius * sin(phi));
        index=index+1;
    end
    signo=signo*-1;
end


