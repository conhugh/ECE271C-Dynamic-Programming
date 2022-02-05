function [moves,flatBoards,boards,stateMap] = tetrisBuild(nr,nc,Pieces,numRots)

%
% Build the basic elements of Tetris game
%
% moves{i}{:} := graphical representation of moves of piece i
% flatboards := row-wise binary representation of boards
% boards{:} := graphical representation of boards
% stateMap := enumerate of states (board # X piece #)
% nr := number of rows of board
% nc := number of columns of board
% Pieces{:} := graphical representation of each piece
% numRots := number of allowable rotations per piece

nPieces = length(Pieces);

disp('Building moves...')

for pNo = 1:nPieces,
    piece = Pieces{pNo};
    kc = 0;
    for rNo = 0:numRots(pNo),
        if rNo > 0,
            piece = rot90(piece);
        end
        [pRow,pCol] = size(piece);
        for k = 1:(nc-pCol+1),
            kc = kc+1;
            blank = zeros(pRow,nc);
            blank(1:pRow,k:(k+pCol-1)) = piece;
            moves{pNo}{kc} = blank;
        end
    end
end

if nargout > 1,
    
    disp('Building boards...')
    
    flatBoards = [];
    nBits = nc*nr;
    nCombos = 2^nBits;
    for kc = 0:(nCombos-1),
        newRow = dec2bin(kc,nBits);
        flatBoards = [flatBoards;newRow];
        boards{kc+1} = reshape(newRow,nr,nc);
    end
    
    stateMap = [];
    for i = 1:nPieces,
        stateMap = [stateMap;[1:nCombos]' i*ones(nCombos,1)];
    end
    
end