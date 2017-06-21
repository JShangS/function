function u = unitfromangles(ang)
%input  ang, an array of n angles
%output u, a unitary n+1 dimensional vector with direction ang

%x = [sin(a)*sin(b)*sin(c)*sin(d) 
%     sin(a)*sin(b)*sin(c)*cos(d) 
%     sin(a)*sin(b)*cos(c)        
%     sin(a)*cos(b)               
%     cos(a) ]                    


N=max(size(ang));
u=[];
for i=1:N    
    u(i)=1;   
    for j=1:N-i+1        
        u(i)=u(i)*sin(ang(j));         
    end
    if (i>1)
        u(i)=u(i)*cos(ang(N-(i-2)));
    end    
end
u(N+1)= cos(ang(1));
end