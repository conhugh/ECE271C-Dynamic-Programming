function [newBoard,score] = nextBoard(board,move)

%
% Computes the next board given current board and move
% Returns empyt matrix if new height exceed nr
%

[nr,nc] = size(board);

top = boardHeight(board);
bot = boardHeight(flip(move));
tmp = find(bot == 0);
bot(tmp) = -inf;

newTop = max(top+bot);

if newTop >= nr,
    bigMove = [move;zeros(newTop-size(move,1),nc)];
    bigBoard = [zeros(newTop-nr,nc);board];
    newBoard = bigMove + bigBoard;
else
    newBoard = board;
    tmp = nr - newTop + 1;
    newBoard(tmp:tmp+size(move,1)-1,:) = ...
        board(tmp:tmp+size(move,1)-1,:) + move;
end

rows = find(sum(newBoard,2) == nc);
score = length(rows);
newBoard(rows,:) = [];

if size(newBoard,1) < nr,
    deficit = nr - size(newBoard,1);
    if deficit > 0,
        newBoard = [zeros(deficit,nc);newBoard];
    end
end
        
