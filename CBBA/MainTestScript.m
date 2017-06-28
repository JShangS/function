% Copyright 2010
% Massachusetts Institute of Technology
% All rights reserved
% Developed by the Aerospace Controls Lab, MIT

%---------------------------------------------------------------------%
% Consensus-Based Bundle Algorithm (CBBA)

% This software package implements the Consensus-Based Bundle Algorithm
% (CBBA), a decentralized market-based protocol that provides provably good 
% approximate solutions for multi-agent multi-task allocation problems
% over networks of heterogeneous agents. The current version supports 
% tasks with time windows of validity, heterogeneous agent-task 
% compatibility requirements, and score functions that balance task
% reward and fuel costs.
%---------------------------------------------------------------------%

%---------------------------------------------------------------------%
% Main test file.  Initializes problem and calls CBBA.
%---------------------------------------------------------------------%

% Clear environment
close all; clear all;
addpath(genpath(cd));

% profile on

SEED = 24377;
rand('seed', SEED);

%---------------------------------------------------------------------%
% Initialize global variables
%---------------------------------------------------------------------%

WORLD.CLR  = rand(100,3);

WORLD.XMIN = -2.0;
WORLD.XMAX =  2.5;
WORLD.YMIN = -1.5;
WORLD.YMAX =  5.5;
WORLD.ZMIN =  0.0;
WORLD.ZMAX =  2.0;
WORLD.MAX_DISTANCE = sqrt((WORLD.XMAX - WORLD.XMIN)^2 + ...
                          (WORLD.YMAX - WORLD.YMIN)^2 + ...
                          (WORLD.ZMAX - WORLD.ZMIN)^2);
 
%---------------------------------------------------------------------%
% Define agents and tasks
%---------------------------------------------------------------------%
% Grab agent and task types from CBBA Parameter definitions
CBBA_Params = CBBA_Init(0,0);

% Initialize possible agent fields
agent_default.id    = 0;            % agent id
agent_default.type  = 0;            % agent type
agent_default.avail = 0;            % agent availability (expected time in sec)
agent_default.clr = [];             % for plotting

agent_default.x       = 0;          % agent position (meters)
agent_default.y       = 0;          % agent position (meters)
agent_default.z       = 0;          % agent position (meters)
agent_default.nom_vel = 0;          % agent cruise velocity (m/s)
agent_default.fuel    = 0;          % agent fuel penalty (per meter)

% FOR USER TO DO:  Set agent fields for specialized agents, for example:
% agent_default.util = 0;

% Initialize possible task fields
task_default.id       = 0;          % task id
task_default.type     = 0;          % task type
task_default.value    = 0;          % task reward
task_default.start    = 0;          % task start time (sec)
task_default.end      = 0;          % task expiry time (sec)
task_default.duration = 0;          % task default duration (sec)
task_default.lambda   = 0.1;        % task exponential discount

task_default.x        = 0;          % task position (meters)
task_default.y        = 0;          % task position (meters)
task_default.z        = 0;          % task position (meters)

% FOR USER TO DO:  Set task fields for specialized tasks

%---------------------------%

% Create some default agents

% QUAD
agent_quad          = agent_default;
agent_quad.type     = CBBA_Params.AGENT_TYPES.QUAD; % agent type
agent_quad.nom_vel  = 2;         % agent cruise velocity (m/s)
agent_quad.fuel     = 1;         % agent fuel penalty (per meter)

% CAR
agent_car           = agent_default;
agent_car.type      = CBBA_Params.AGENT_TYPES.CAR;  % agent type
agent_car.nom_vel   = 2;         % agent cruise velocity (m/s)
agent_car.fuel      = 1;         % agent fuel penalty (per meter)

% Create some default tasks

% Track
task_track          = task_default;
task_track.type     = CBBA_Params.TASK_TYPES.TRACK;      % task type
task_track.value    = 100;    % task reward
task_track.start    = 0;      % task start time (sec)
task_track.end      = 100;    % task expiry time (sec)
task_track.duration = 5;      % task default duration (sec)

% Rescue
task_rescue          = task_default;
task_rescue.type     = CBBA_Params.TASK_TYPES.RESCUE;      % task type
task_rescue.value    = 100;    % task reward
task_rescue.start    = 0;      % task start time (sec)
task_rescue.end      = 100;    % task expiry time (sec)
task_rescue.duration = 15;     % task default duration (sec)


%---------------------------------------------------------------------%
% Define sample scenario
%---------------------------------------------------------------------%

N = 5;      % # of agents
M = 10;     % # of tasks

% Create random agents
for n=1:N

    if(n/N <= 1/2)
        agents(n) = agent_quad;
    else
        agents(n) = agent_car;
    end

    % Init remaining agent params
    agents(n).id   = n;
    agents(n).x    = rand(1)*(WORLD.XMAX - WORLD.XMIN) + WORLD.XMIN;
    agents(n).y    = rand(1)*(WORLD.YMAX - WORLD.YMIN) + WORLD.YMIN;
    agents(n).clr  = WORLD.CLR(n,:);
end

% Create random tasks
for m=1:M
    if(m/M <= 1/2)
        tasks(m) = task_track;
    else
        tasks(m) = task_rescue;
    end
    tasks(m).id       = m;
    tasks(m).start    = rand(1)*100;
    tasks(m).end      = tasks(m).start + 1*tasks(m).duration;
    tasks(m).x        = rand(1)*(WORLD.XMAX - WORLD.XMIN) + WORLD.XMIN;
    tasks(m).y        = rand(1)*(WORLD.YMAX - WORLD.YMIN) + WORLD.YMIN;
    tasks(m).z        = rand(1)*(WORLD.ZMAX - WORLD.ZMIN) + WORLD.ZMIN;
end

%---------------------------------------------------------------------%
% Initialize communication graph and diameter
%---------------------------------------------------------------------%

% Fully connected graph
Graph = ~eye(N);

%---------------------------------------------------------------------%
% Run CBBA
%---------------------------------------------------------------------%

[CBBA_Assignments Total_Score] = CBBA_Main(agents, tasks, Graph)

PlotAssignments(WORLD, CBBA_Assignments, agents, tasks, 1);


% profile off
% profile report


