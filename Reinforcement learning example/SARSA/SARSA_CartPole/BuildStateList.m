function [ states ] = BuildStateList
%BuildStateList builds a state list from a state matrix

% state discretization for the mountain car problem
x1div = (2-(-2)) / 3.0;
x2div = (0.1-(-0.1)) / 2.0;
x3div = (deg2rad(12)-(deg2rad(-12)))/8.0;
x4div = (deg2rad(10)-(deg2rad(-10)))/2.0;

%x1  = -2:x1div:2;
x1  = [-1  1];
%x2  = -0.5:x2div:0.5;
x2  = [-1 0 1];
x3  = deg2rad(-12):x3div:deg2rad(12);
%x4  = deg2rad(-10):x4div:deg2rad(10);
x4   = [-1  1];


I=size(x1,2);
J=size(x2,2);
K=size(x3,2);
L=size(x4,2);
states=[];
index=1;
for i=1:I    
    for j=1:J
        for k=1:K
            for l=1:L
                states(index,1)=x1(i);
                states(index,2)=x2(j);
                states(index,3)=x3(k);
                states(index,4)=x4(l);
                index=index+1;
            end
        end
    end
end
