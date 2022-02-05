% Homework 2 Matlab Problems 5 - Deterministic
% Max Emerick

% Start by running setup from tetrisMain

clc; clear; close all

tData.Deterministic = 1; % is the next piece random or determinist
                         % random -> 0
                         % deterministic start with piece # -> # in {1,2,3}

tData.nMaxPieces = 100; % max number of pieces per episode
tData.nEpisodes = 1; % number of episodes

tData.buildStates = 1; % flag to build state space

tData.GameSize = [6,3]; % height x width
tData.RowCap = 3; % height of gameOver
tData.TimeDelay=.05; % Time delay of dropping piece (lower number=faster)

tData.SelectMethod = 1; %1. Value iteration.  2. Q-Learning Exact

% Define the game pieces
Pieces{1} = [0 1;1 1];
Pieces{2} = [0 1 1;1 1 0];
Pieces{3} = [1 ; 1];

numRots = [3 1 1]; %number of times each piece can be rotated
    
tData.Pieces = Pieces;
tData.numRots = numRots;
% Set piece colors
tData.Pcolor=1:length(Pieces);

% Build tetris states & moves
if tData.buildStates
    display(datetime('now'))
    tic
    [moves,flatBoards,boards,stateMap] = ...
        tetrisBuild(tData.RowCap,tData.GameSize(2),Pieces,numRots);
    tData.flatBoards = flatBoards;
    tData.boards = boards;
    tData.moves = moves;
    tData.stateMap = stateMap;
    toc
else
    moves = tetrisBuild(tData.RowCap,tData.GameSize(2),Pieces,numRots);
    tData.moves = moves;
end

tData.S_Sounds=0; % Switch, 1=sounds on, 0=sounds off
tData.S_Plot=1;  % Switch to Perform plotting, 1=yes

if tData.Deterministic == 0
    tData.startPiece = randi(length(Pieces));
else
    tData.startPiece = tData.Deterministic;
end

% Don't Run Demo
%[Iscore,nPieces] = tetrisMyPlayDemo(tData);

% Initialize cost matrix
numStates = length(stateMap(:,1));
numBoards = length(boards);
M = Inf*ones(numStates,numBoards+1); % stores costs
C = zeros(numStates,numBoards+1); % stores controls

% Populate cost matrix
for i = 1:numStates % for each state
    
    board = boards{mod(i-1,numBoards)+1}; % get associated board
    pieceNum = stateMap(i,2); % get associated piece number
    
    for u = 1:length(tData.moves{pieceNum}) % for each admissible control
        move = tData.moves{pieceNum}{u}; % get move
        [newBoard,score] = nextBoard(board,move); % find next state and score
        height = max(boardHeight(newBoard)); % check if game ended
        if height > tData.RowCap
            j = numBoards+1; % if game ended, transition to end game state
        else % index of next board
            j = b2d(flip(reshape(newBoard,1,numel(newBoard))))+1;
        end
        M(i,j) = -score; % cost = negative of rows eliminated
        C(i,j) = u;
    end

end

% Setup information
piece1 = 3;
numStages = 100;

% Initialize cost tensor
L = Inf*ones(numBoards+1,numBoards+1,numStages);

for i = 1:numStages % for each stage
    
    % Set stage cost matrix according to correct subset of state space
    piece = mod(i+piece1+1,3)+1;
    base_index = numBoards*(piece-1);
    L(1:numBoards,:,i) = M(base_index+1:base_index+numBoards,:);
    
end
L(numBoards+1,numBoards+1,:) = zeros(1,numStages); % no cost after game ends


% Solve shortest path problem
[a,b,c] = size(L); % store size of cost tensor
cost_to_go = Inf*ones(a,c+1); % initialize array to store cost-to-go for each state
cost_to_go(:,c+1) = zeros(a,1); % set cost-to-go for terminal state to zero
next_node = zeros(a,c); % array holding optimal transition for each state
for k = c:-1:1 % for each stage
    for i = 1:a % for each state in stage
        total_cost_i = Inf*ones(1,b); % initialize vector of costs for each control
        for j = 1:b % for each state in next stage
            stage_cost = L(i,j,k); % get stage cost-to-go for next state
            total_cost_i(j) = stage_cost + cost_to_go(j,k+1); % optimal cost through node j
        end
        [cost,node] = min(total_cost_i); % find optimal cost-to-go for node i
        cost_to_go(i,k) = cost; % store cost-to-go
        next_node(i,k) = node; % store corresponding optimal control
    end
end

route = zeros([1,k]); % initialize vector to store optimal path
%[cost,node1] = min(cost_to_go(:,1)); % find optimal first step

% Enforce start at empty board
route(1) = 1;
cost = cost_to_go(1,1);
control = zeros([1,k-1]);

for k = 2:c+1 % for each remaining stage
    node = next_node(route(k-1),k-1); % follow optimal transition to next state
    route(k) = node; % add optimal next state to path
end

% Construct and save optimal control sequence
for i = 1:c
    
    piece = mod(i+piece1+1,3)+1;
    base_index = numBoards*(piece-1);

    control(i) = C(route(i)+base_index,route(i+1));
    
end

save control


