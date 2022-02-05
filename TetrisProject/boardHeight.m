function h = boardHeight(board)

%
% Compute the height of each column
%

[nr,nc] = size(board);

for i = 1:nc,
    fboard = flip(board);
    tmp = max(find(fboard(:,i) == 1));
    if isempty(tmp),
        top(i) = 0;
    else
        top(i) = tmp;
    end
end

h = top;