%Connor Hughes
%ECE 271C HW3

%% Problem 6: Backgammon
n = 6;
k = 8;
numDice = 2;

%generate possible initial board configurations for n spaces, k pieces
combs = nmultichoosek(linspace(n, 1, n), k);
%generate the possible unique (order-indep.) rolls of dice
rolls = nmultichoosek([1, 2, 3, 4, 5, 6], numDice);
%generate cell structure to store state information
states = cell(length(combs(:, 1))*length(rolls(:, 1)), 2);

%generate states:
for i = 1:length(states)
   %reformat each board config as a 6-element array,
   %where each element gives piece count at corresponding spot
   [GC, GR] = groupcounts(combs(floor((i - 1)/21) + 1, :)');
   states{i, 1} = zeros(1, n); 
   for j = 1:length(GR)
      states{i, 1}(1, GR(j)) = GC(j); 
   end
   states{i, 2} = rolls(mod(i - 1, 21) + 1, :);
end

%generate control set for a given state:
%test = max(find(states{27022, 1}))
% test0 = find(states{4066, 1}(1, 3:end))
% test = length(find(states{4066, 1}(1, 3:end)))
test6 = getSingleMoves(states{4066, 1}, 6);
test5 = getSingleMoves(states{4066, 1}, 5);
test4 = getSingleMoves(states{4066, 1}, 4);
test3 = getSingleMoves(states{4066, 1}, 3);
test2 = getSingleMoves(states{4066, 1}, 2);
test1 = getSingleMoves(states{4066, 1}, 1);

%% Problem 7 (Tetris) Setup
% Define the game pieces
Pieces{1} = [0 1;1 1];
Pieces{2} = [0 1 1;1 1 0];
Pieces{3} = [1 ; 1];
numRots = [3 1 1]; %number of times each piece can be rotated
%P = [0 1 0; 0 0 1; 1 0 0]; %piece transition probability matrix
P = [1/3 1/3 1/3; 1/3 1/3 1/3; 1/3 1/3 1/3];

[moves, flatBoards, boards, stateMap] = tetrisBuild(3, 3, Pieces, numRots);

%generate network of states for 100 stages
stageMap = cell(1, 100);
for k = 1:100
   stageMap{k} = stateMap;
end

gameover = [NaN, 1];
stateMap = [stateMap; gameover];

%test nextBoard.m:
[newBoard, score] = nextBoard(boards{1}, moves{1}{2});
[newerBoard, score2] = nextBoard(newBoard, moves{3}{1});

%generate sequence of pieces (deterministic case):
pStartNum = 3;
pieces = zeros(1, 101);
pieces(1) = pStartNum;
for k = 2:100
   pNum = nextPiece(P, pieces(k - 1));
   pieces(k) = pNum;
end
pieces(101) = 1;



%% Problem 7 (Tetris) CTG Matrix Generation
for totalPcs = 1:200
    %generate sequence of pieces (deterministic case):
    pStartNum = 3;
    pieces = zeros(1, 101);
    pieces(1) = pStartNum;
    for k = 2:100
       pNum = nextPiece(P, pieces(k - 1));
       pieces(k) = pNum;
    end
    pieces(101) = 1;

    %generate single stage reward matrix:
    rtg = cell(1537, 101);
    for k = 1:100
        for i = 1:1537
            rtg{i, k} = NaN(1, 4);
        end
    end
    for i = 1:1537
        rtg{i, 101} = zeros(1, 4);
    end
     for k = 100:-1:1
        %for each (non-gameover) state at current stage
        for i = 1:1536
            %get current board:
            currBoard = boards{mod((i - 1), 512) + 1};
            %initialize options matrix:
            opts = NaN(length(moves{stateMap(i, 2)}), 3);
            %for each possible move with the piece assoc'd with current state
            for j = 1:length(moves{stateMap(i, 2)})
                move = moves{stateMap(i, 2)}{j};
                opts(j, 1) = j;
                [newBoard, score] = nextBoard(currBoard, move);  %generate next board
                if length(newBoard(:, 1)) <= 3
                    %find states at next stage which have the right next board
                    boardIndices = [board2dec(newBoard) + 1, board2dec(newBoard) + 513, board2dec(newBoard) + 1025];
                    %find associated rewards-to-go for each next state option
                    rtgOpts = [rtg{boardIndices(1), k + 1}(2), rtg{boardIndices(2), k + 1}(2), rtg{boardIndices(3), k + 1}(2)];
                    %compute expected value for reward-to-go for next state:
                    expRtg = P(stateMap(i, 2), 1)*rtgOpts(1) + P(stateMap(i, 2), 2)*rtgOpts(2) + P(stateMap(i, 2), 3)*rtgOpts(3);
                    %boardInd = boardIndices(nextPiece(P, stateMap(i, 2)));
                    boardInd = boardIndices(pieces(k + 1));
                else
                    boardInd = 1537;  %set index of next board to gameover
                    score = 0; %store reward for getting to next board
                    expRtg = 0;
                end
                opts(j, 2) = score + expRtg; %compute expected rtg 
                opts(j, 3) = boardInd;  %store next board index
                opts(j, 4) = stateMap(boardInd, 2); %store piece for next state
            end
            [maxRwd, moveNum] = max(opts(:, 2));  %select max expected rtg option
            rtg{i, k} = [moveNum, maxRwd, opts(moveNum, 3), pieces(k + 1)]; %store move, rtg, next board index, next piece
        end
        rtg{1537, k} = [NaN, 0, 1537, 1]; %set gameover state for this stage
     end
 
end

%% Problem 7 (Tetris) Optimal Path Assembly
%start with empty board, and 1 for pc 1, 513 for pc 2, 1025 for pc 3
currState = 1;
path = zeros(1, 100);
boardTracker = zeros(1, 100);
%follow sequence of next states through subsequent stages
for k = 1:100
   %store the move for each state
   boardTracker(k) = currState;
   path(k) = rtg{currState, k}(1);
   nextState = rtg{currState, k}(3);
   currState = nextState;
end
save path

%% Helper Functions
%n choose k with replacement
function combs = nmultichoosek(values, k)
    %// Return number of multisubsets or actual multisubsets.
    if numel(values)==1
        n = values;
        combs = nchoosek(n+k-1,k);
    else
        n = numel(values);
        combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
        combs = reshape(values(combs),[],k);
    end
end
%generate set of legal moves for single die roll
function nextBoards = getSingleMoves(board, die)
    %find max board position with nonzero value
    maxPos = max(find(board));
    if die >= maxPos
        %can only move highest-position piece to home
        %generate possible next board
        diff = zeros(1, 6);
        diff(max(find(board))) = -1;
        nextBoards = board + diff;
    else 
        %can move any piece as long as it doesn't go to home
        next = zeros(length(find(board(1,(die + 1):end))), 6);
        ind = 1;
        %generate possible next boards
        for pos = (die + 1):max(find(board))
           if board(pos) ~=0
               diff = zeros(1, 6);
               diff(pos) = -1;
               diff(pos - die) = 1;
               next(ind, :) = board + diff;
               ind = ind + 1;
           end
        end
        nextBoards = next;
    end  
end
%generate set of legal moves for roll of two dice
function [nextBoards, moves] = getDoubleMoves(board, dice)
    
end
%random piece selection given current piece and transition prob matrix
function pc = nextPiece(P, currPc)
    CDF = zeros(1, 3);
    CDF(1) = P(currPc, 1);
    CDF(2) = P(currPc, 1) + P(currPc, 2);
    CDF(3) = P(currPc, 1) + P(currPc, 2) + P(currPc, 3);
    n = rand;
    pc = 0;
    if n <= CDF(1)
        pc = 1;
    else
        if n <= CDF(2)
            pc = 2;
        else
            pc = 3;
        end
    end
end

%convert board to decimal value
function[dec] = board2dec(board)
    decSum = 0;
    for i = 1:length(board(:, 1))
        for j = 1:length(board(1, :))
            if board(i, j) == 1
                decSum = decSum + 2^(3*(j - 1) + i - 1);
            end
        end
    end
    dec = decSum;
end