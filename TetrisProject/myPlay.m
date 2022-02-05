function [decision,DATAout] = myPlay(board,pieceNum,DATA)
% Testing out various policies
% 
% Places the piece to minimize the maximum height
%

height = NaN(1,length(DATA.moves{pieceNum}));

squares = NaN(1,length(DATA.moves{pieceNum}));

holes = NaN(1,length(DATA.moves{pieceNum}));

avg_diffs = NaN(1,length(DATA.moves{pieceNum}));

scores = NaN(1,length(DATA.moves{pieceNum}));

for u = 1:length(DATA.moves{pieceNum})
    move = DATA.moves{pieceNum}{u};
    [theNextBoard,nextScore] = nextBoard(board,move);
    
    heights = boardHeight(theNextBoard);
    
    diffs = heights(2:length(heights))-heights(1:length(heights)-1);
    avg_diffs(u) = mean(abs(diffs));
    
    height(u) = max(heights);
    squares(u) = sum(heights.^2);
    
    scores(u) = nextScore;
    
    holes(u) = sum(heights - sum(theNextBoard,1));
    
end

% Eliminate row if possible, place low if not
% [maxScore,uBest] = max(nextScore);
% if maxScore == 0
%     [~,uBest] = min(height);
% end

% minimize height
%[~,uBest] = min(height);

% minimize square of column heights
%[~,uBest] = min(squares);

% Minimize trapped squares, then square of heights
%[~,uBest] = min(squares+1000*holes);

% Minimize jaggedness, then square of heights
%[~,uBest] = min(squares+1000*avg_diffs);

% Weighted Average - not best
%[~,uBest] = min(squares+0.2*holes);

% Random
uBest = randi(length(height));

decision = DATA.moves{pieceNum}{uBest};
DATAout = DATA;
