function [ states ] = BuildStateList
%BuildStateList builds a state list from a state matrix

x1 = [-0.5 0 0.5]; %difference in x
x2 = [-0.5 0 0.5]; %difference in y
%x3 = [-1 0 1]; % x coordinate



I=size(x1,2);
J=size(x2,2);
%K=size(x3,2);
%L=size(x4,2);
states=[];
index=1;
for i=1:I    
    for j=1:J
        %for k=1:K           
            %for l=1:L
                    states(index,1)=x1(i);
                    states(index,2)=x2(j);
                    %states(index,3)=x3(k);               
                    %states(index,4)=x4(l);               
                    index=index+1;  
            %end         
        %end        
    end
end
