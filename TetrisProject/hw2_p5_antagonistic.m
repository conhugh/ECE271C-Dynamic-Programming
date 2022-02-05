% Homework 2 Matlab Problems 5 - Antagonistic
% Max Emerick

% Start by running setup from tetrisMain

clc; clear; close all;

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

% Setup information
piece1 = 1;
numStages = 100;
numStates = length(stateMap(:,1));
numBoards = length(boards);


% Find costs and solve shortest path problem simultaneously
a = numStates+1;
cost_to_go = Inf*ones(a,numStages+1); % initialize array to store cost-to-go for each state
cost_to_go(:,numStages+1) = zeros(a,1); % set cost-to-go for terminal state to zero
next_node = zeros(a,numStages); % array holding expected transition for each state
best_control = zeros(a,numStages); % array holding best control for each state

for k = numStages:-1:1 % for each stage
    for i = 1:a % for each state in stage
        
        if i == 1537 % if game is over, set values acoordingly
            cost_to_go(i,k) = cost_to_go(i,k+1);
            best_control(i,k) = 0;
            next_node(i,k) = 1537;
        else % game is not over
            board = boards{mod(i-1,numBoards)+1}; % get associated board
            pieceNum = stateMap(i,2); % get associated piece number
            numMoves = length(tData.moves{pieceNum}); % find number of controls
            
            % Initialize matrices to store min-max costs and expected transition
            mmatrix = Inf*ones(3,numMoves); % stores costs
            expected_node = zeros(3,numMoves); % stores expected transitions
            max_costs = Inf*ones(1,numMoves); % stores worst-case cost for each control

            for u = 1:numMoves % for each admissible control
                move = tData.moves{pieceNum}{u}; % get move
                [newBoard,score] = nextBoard(board,move); % find next state and score
                height = max(boardHeight(newBoard)); % check if game ended
                jx = 1537; % initialize next-node variables to game over
                j = 1537; % figure out what's going on here
                if height > tData.RowCap
                    jx = numStates+1; % if game ended, transition to end game state
                else % get index of next board arrangment
                    jx = b2d(flip(reshape(newBoard,1,numel(newBoard))))+1;
                end

                for x = 1:3 % for each possible next piece
                    if jx ~= numStates+1 % check if already in game over
                        j = jx+(x-1)*numBoards; % if not find next state
                    end
                    
                    % find cost-to-go through next state given next piece
                    mmatrix(x,u) = -score + cost_to_go(j,k+1);
                    expected_node(x,u) = j; % store expected transition

                end
                
                % find and store worst-case cost for each control
                max_costs(u) = max(mmatrix(:,u));

            end
            
            % find and store cost, control, and next state corresponding to least worst-case
            [mcost,control] = min(max_costs);
            [asdf,next_piece] = max(mmatrix(:,control));
            cost_to_go(i,k) = mcost;
            best_control(i,k) = control;
            next_node(i,k) = expected_node(next_piece,control);
            
        end
        
    end
    
end

% initialize vector to store optimal path
route = zeros([1,numStages+1]);

% Enforce start at empty board given first piece
route(1) = numBoards*(piece1-1)+1;
cost = cost_to_go(route(1),1);
control = zeros([1,numStages]);

for k = 2:numStages+1 % for each remaining stage
    node = next_node(route(k-1),k-1); % follow optimal transition to next state
    route(k) = node; % add optimal next state to path
end

% Save optimal policy for simulation
save best_control;


