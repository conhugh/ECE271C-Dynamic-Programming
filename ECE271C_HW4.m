%Connor Hughes
%ECE 271C HW4

%% Problem 2
%initialize parameter struct:
p.t = 0.2;    %likelihood of transition good->bad, assuming starts good
p.I = 1;      %cost of inspection
p.R = 3;      %cost of replacement
p.C = 2;      %cost of producing bad item
p.N = 8;

x0 = [1, 0];  %initial belief - machine is in good state

%Initialize cost-to-go matrix, next node matrix, states matrix:
S = cell(1, p.N + 1);
CTG = cell(1, p.N + 1);
Paths = cell(1, p.N + 1);
NN = cell(1, p.N + 1);
for k = 1:(p.N + 1)
    CTG{k}= NaN(k, 1);
    Paths{k}= NaN(k, p.N + 2 - k);
    S{k}= NaN(k, 2);
    NN{k} = NaN(k, 1);
end
CTG{p.N + 1}(:) = 0;
NN{p.N + 1}(:) = NaN;

%populate states matrix:
S{1}(1, :) = x0;
for k = 2:(p.N + 1) 
    S{k}(1, :) = [1, 0];
    for i = 2:length(S{k})
        S{k}(i, :) = [0.8*S{k - 1}(i - 1, 1), 0.2*S{k - 1}(i - 1, 1) + S{k - 1}(i - 1, 2)];
    end
end

%populate cost-to-go matrix and paths matrix:
for k = (p.N):-1:1
    for i = 1:length(CTG{k})
        opts = NaN(2, 1);
        opts(1, 1) = gk(S{k}(i, :), 1, p) + CTG{k + 1}(1);
        opts(2, 1) = gk(S{k}(i, :), 0, p) + CTG{k + 1}(i + 1);
        [m, ind] = min(opts);
        CTG{k}(i) = m;
        if ind == 1
            NN{k}(i) = 1;
        end
        if ind == 2
            NN{k}(i) = i + 1;
        end
    end
end

%% Problem 4
%see the last section of this script for function implementation
S = inventory(5, 3, 3, 0.5, 2, 1, 1.5, 5);

%% Problem 5
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

maxMvs = 200;
%generate sequence of pieces (deterministic case):
pStartNum = 3;
pieces = zeros(1, maxMvs + 1);
pieces(1) = pStartNum;
for k = 2:maxMvs
   pNum = nextPiece(P, pieces(k - 1));
   pieces(k) = pNum;
end
pieces(maxMvs + 1) = 1;

%generate single stage reward matrix:
rtg = cell(1537, maxMvs + 1);
for k = 1:maxMvs
    for i = 1:1537
        rtg{i, k} = NaN(1, 4);
    end
end
for i = 1:1537
    rtg{i, maxMvs + 1} = zeros(1, 4);
end
 for k = maxMvs:-1:1
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

Jstar = NaN(1, 200);
for i = 200:-1:1
   Jstar(201 - i) = rtg{1025, i}(2);
end
save('RTG', 'rtg')

%% Problem 5 plots
close all
load JstarPc1.mat
stg = linspace(1, 200, 200);
scatter(stg, Jstar, 'Marker', '.', 'LineWidth', 2)
xlabel("Number of Stages")
ylabel("Cost to Go")
hold on
load JstarPc2.mat
scatter(stg, Jstar, 'Marker', '.', 'LineWidth', 2)
hold on
load JstarPc3.mat
scatter(stg, Jstar, 'Marker', '.', 'LineWidth', 2)
legend('Piece 1', 'Piece 2', 'Piece 3')
title('Expected No. Rows Eliminated vs. Number of Stages Allowed')
ax = gca
ax.FontSize = 18
%% Helper Functions and Problem 4
%Problem 2 stage cost:
function cst = gk(xk, uk, p)
if uk == 1
   cst = p.I*xk(1) + (p.C + p.I + p.R)*xk(2);
end
if uk == 0
    cst = p.C*xk(2);
end

end

function S = inventory(Bmax, Bmin, wmax, q0, p, h, c, N)
    %initialize network of states
    nStates = Bmax + abs(Bmin) + 1;
    States = zeros(nStates, N);
    for k = 1:N
       States(:, k) = linspace(-Bmin, Bmax, nStates); 
    end

    %initialize CTG matrix
    CTG = NaN(nStates, N);
    CTG(:, N) = zeros(nStates, 1);

    %initialize matrix of next nodes (for assembling optimal path later)
    %NN = NaN(nStates, N);

    %initialize matrix of controls
    U = NaN(nStates, N);
    
    %populate CTG matrix
    for k = (N - 1):-1:1
        for i = 1:nStates
            xk = States(i, k);
            Umax = Bmax - xk + wmax;
            J = Inf(Umax + 1, 1);
            for u = 0:Umax
                if(round(xk + u - (1-q0)*wmax) >= Bmax)
                    J(u + 1) = c*u + p*max(0, -(xk + u - (1-q0)*wmax)) + h*max(0,(xk + u - (1-q0)*wmax)) + CTG(find(States(:, k + 1)==Bmax, 1), k + 1);    
                end
                if (round(xk + u - (1-q0)*wmax) <= -Bmin)
                    J(u + 1) = c*u + p*max(0, -(xk + u - (1-q0)*wmax)) + h*max(0,(xk + u - (1-q0)*wmax)) + CTG(find(States(:, k + 1)==(-Bmin), 1), k + 1);
                else
                    J(u + 1) = c*u + p*max(0, -(xk + u - (1-q0)*wmax)) + h*max(0,(xk + u - (1-q0)*wmax)) + CTG(find(States(:, k + 1)==round(xk + u - (1-q0)*wmax), 1), k + 1);    
                end
            end
            [Jk, ind] = min(J);
            U(i, k) = ind - 1;
            CTG(i, k) = Jk;
            NN(i, k) = find(round((xk + u - (1-q0)*wmax)));
        end
    end
    
    %search the control table to locate thresholds
    T = zeros(N, 1);
    for k = (N - 1):-1:1
        thresh = find(U(:, k) == 0, 1);
        T(k, 1) = States(thresh, k);
    end
    S = T;
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