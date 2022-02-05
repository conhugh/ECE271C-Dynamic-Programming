% tetrisMain.m script
% Enter the values to setup tetris simulation
%Readme%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Definition:
% boardHeight: takes board as an input and return heights of each column
% myPlay: play rule to minimize maximum height, evaluate all possible 
%         control (moves) givrn a piece and returns the one resulting in
%         minimum maximum hight for next board.
% nextBoard: takes board and selected control (move) and returns the new
%            board and positive score if any lines were wiped out.
% randiP: selects a random integer according to probability vector p.
% tetrisBuild: build the basic elements of Tetris game
% tetrisBuildDemo: illustrates the use/structure of different 
%                  functions/variables in a the code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all; 

tData.Deterministic = 0; % is the next piece random or deterministic
                         % random -> 0
                         % deterministic start with piece # -> # in {1,2,3}

tData.nMaxPieces = 1000; % max number of pieces per episode
tData.nEpisodes = 200; % number of episodes

tData.buildStates = 0; % flag to build state space

tData.GameSize = [18,10]; % height x width
tData.RowCap = 15; % height of gameOver
tData.TimeDelay=0; % Time delay of dropping piece (lower number=faster)

tData.SelectMethod = 1; %1. Value iteration.  2. Q-Learning Exact

% Define the game pieces
Pieces{1} = [1 1 0;0 1 1];
Pieces{2} = [0 1 1;1 1 0];
Pieces{3} = [1;1;1;1];
Pieces{4} = [1 1;1 1];
Pieces{5} = [0 1 0;1 1 1];
Pieces{6} = [0 0 1;1 1 1];
Pieces{7} = [1 0 0;1 1 1];

numRots = [1 1 1 0 3 3 3]; %number of times each piece can be rotated
    
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
tData.S_Plot=0;  % Switch to Perform plotting, 1=yes

if tData.Deterministic == 0
    tData.startPiece = randi(length(Pieces));
else
    tData.startPiece = tData.Deterministic;
end

% Initialize parameters for learning
global WEIGHTS;
global WEIGHTUPDATE;
global TIME;
global TRACKWEIGHTS;
TIME = 100; % start larger, change in 2 places
WEIGHTS = [10,10]; % close enough for good play?
WEIGHTUPDATE = WEIGHTS;
TRACKWEIGHTS = NaN(2,100);
    
% Run Demo 
[Iscore,nPieces] = tetrisMyPlayDemo(tData);

%% Data Analysis

%save Iscore;
plot(Iscore);
hold on;
run_avg = cumsum(Iscore)./(1:length(Iscore));
plot(run_avg);
title('Strategy 3 - Running Average');
xlabel('Trial');
ylabel('Score');


% figure;
% histogram(Iscore);
% title('Strategy 3 - Score Distribution');
% xlabel('Score');
% ylabel('Frequency');

figure;
plot(TRACKWEIGHTS(1,:));
hold on;
plot(TRACKWEIGHTS(2,:));
% plot(TRACKWEIGHTS(3,:));
% plot(TRACKWEIGHTS(4,:));
% plot(TRACKWEIGHTS(5,:));
legend('1','2','3','4','5');
legend('1','2');

mean(Iscore)