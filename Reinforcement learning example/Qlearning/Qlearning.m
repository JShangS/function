function robot=Qlearning(robot,task,Rew,N)
% Q learning Arithmatic
[m,n]=size(Rew);           % load the task and the robot
gam=robot.gamma;
alp=robot.alpha;
as=1:4;
si=task.initialState;ki=sub2ind([m,n],si(1),si(2));
st=task.terminalState;kt=sub2ind([m,n],st(1),st(2));

for lp=1:N
    alp=alp*exp((1-rem(lp,10))/N);
    Q=robot.Qtable;        % load the new table
    robot.state=si;        % return to the initial state 
    s0=si;k0=ki;
    step=0;
    slist=[];
    while k0~=kt 
        qs=Q(k0,:);
        q=min(qs);ras=as;ras(qs==q)=[];
        if isempty(ras)    % policy
            a=unidrnd(4,1);
        else
            qs(qs==q)=[]; 
            a=drnd(ras,qs-q(1));               
        end
        s=exact(s0,a); 
        if any(s<1)||any(s>m)||any(s>n)    % beyond the environment
            Q(k0,a)=Q(k0,a)-100; break;
        else               % feasible
            k=sub2ind([m,n],s(1),s(2));                       
            v=max(Q(k,:));                 % V(s) 
            D=Rew(k)+gam*v-Q(k0,a)-.2;
            Q(k0,a)=Q(k0,a)+alp*D;
        end       
        slist=[slist;s0,a,s];
        s0=s;k0=k;step=step+1;
%        disp([num2str(s0),' | ',num2str(a),' | ',num2str(s)])  % ´òÓ¡×´Ì¬×ªÒÆ
    end
    
    robot.state=s;         % update the robot
    robot.Qtable=Q;
    if k0==kt
        if isempty(robot.best), % record
            robot.best=slist;
        elseif step<size(robot.best,1),
            robot.best=slist;
        end  
    end     
end

function s=exact(s,a)
switch(a)
    case 1   % dowm
        s=s+[1,0];
    case 2   % right
        s=s+[0,1];
    case 3   % up
        s=s+[-1,0];
    case 4   % left
        s=s+[0,-1];
    otherwise
end
