% Copyright 2010
% Massachusetts Institute of Technology
% All rights reserved
% Developed by the Aerospace Controls Lab, MIT

%---------------------------------------------------------------------%
% Main CBBA Function
%---------------------------------------------------------------------%

function [CBBA_Data Total_Score] = CBBA_Main(agents, tasks, Graph)

% Initialize CBBA parameters
CBBA_Params = CBBA_Init(length(agents),length(tasks));

for n=1:CBBA_Params.N,
    CBBA_Data(n).agentID          = agents(n).id;
    CBBA_Data(n).agentIndex       = n;
    CBBA_Data(n).bundle           = -ones(1, CBBA_Params.MAX_DEPTH);
    CBBA_Data(n).path             = -ones(1, CBBA_Params.MAX_DEPTH);
    CBBA_Data(n).times            = -ones(1, CBBA_Params.MAX_DEPTH);
    CBBA_Data(n).scores           = -ones(1, CBBA_Params.MAX_DEPTH);
    CBBA_Data(n).bids             = zeros(1, CBBA_Params.M);
    CBBA_Data(n).winners          = zeros(1, CBBA_Params.M);
    CBBA_Data(n).winnerBids       = zeros(1, CBBA_Params.M);
end

% Initialize working variables
T         = 1;                                      % Current iteration
t         = zeros(CBBA_Params.N, CBBA_Params.N);    % Matrix of time of updates from the current winners
lastTime  = T-1;
doneFlag  = 0;

% Main CBBA loop (runs until convergence)
while(doneFlag == 0)
    
    %---------------------------------------%
    % 1. Communicate
    %---------------------------------------%
    % Perform consensus on winning agents and bid values (synchronous)
    [CBBA_Data t] = CBBA_Communicate(CBBA_Params, CBBA_Data, Graph, t, T);

    %---------------------------------------%
    % 2. Run CBBA bundle building/updating
    %---------------------------------------%
    % Run CBBA on each agent (decentralized but synchronous)
    for n=1:CBBA_Params.N,

        [CBBA_Data(n) newBid] = CBBA_Bundle(CBBA_Params, CBBA_Data(n), agents(n), tasks);

        % Update last time things changed (needed for convergence but will
        % be removed in the final implementation)
        if(newBid),
            lastTime = T;
        end
    end

    %---------------------------------------%
    % 3. Convergence Check
    %---------------------------------------%
    % Determine if the assignment is over (implemented for now, but later
    % this loop will just run forever)
    if(T-lastTime > CBBA_Params.N)
        doneFlag   = 1;
    elseif(T-lastTime > 2*CBBA_Params.N)
        disp('Algorithm did not converge due to communication trouble');
        doneFlag = 1;
    else
        % Maintain loop
        T = T + 1;
    end
end


% Map path and bundle values to actual task indices
for n=1:CBBA_Params.N,
    for m=1:CBBA_Params.MAX_DEPTH,
        if(CBBA_Data(n).bundle(m) == -1),
            break;
        else
            CBBA_Data(n).bundle(m) = tasks(CBBA_Data(n).bundle(m)).id;
        end
        
        if(CBBA_Data(n).path(m) == -1),
            break;
        else
            CBBA_Data(n).path(m) = tasks(CBBA_Data(n).path(m)).id;
        end
    end
end

% Compute the total score of the CBBA assignment
Total_Score = 0;
 for n=1:CBBA_Params.N,
     for m=1:CBBA_Params.MAX_DEPTH,
         if CBBA_Data(n).scores(m) > -1,
             Total_Score = Total_Score + CBBA_Data(n).scores(m);
         else
             break;
         end
     end
 end
