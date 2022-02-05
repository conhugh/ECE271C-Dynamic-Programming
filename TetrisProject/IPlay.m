function [decision,DATAout] = IPlay(board,pieceNum,DATA)
% 
% Places the piece where you want 
% Use left/right arrow to cycle through actions and down arrow to place

Moves = length(DATA.moves{pieceNum});
go = 0;

current = 1;

while go == 0
    move = DATA.moves{pieceNum}{current};
    [r,c] = size(board);  % Get the matrix size
    Newboard = [move;zeros(3-size(move,1),c);board];
    [r,c] = size(Newboard);  % Get the matrix size
    imagesc((1:c)+0.5,(1:r)+0.5,Newboard,[0,max(DATA.Pcolor)]);            % Plot the image
    axis equal                                   % Make axes grid sizes equal
    set(gca,'XTick',1:c+1,'YTick',1:r+1,...  % Change some axes properties
    'XLim',[1 c+1],'YLim',[1 r+1],...
    'XTickLabel',0:c,'YTickLabel',0:r,...
    'GridLineStyle','-','XGrid','on','YGrid','on');
    hold on
    plot([0,r],[r-DATA.RowCap+1,r-DATA.RowCap+1],'r', 'LineWidth',3)
    hold off

    w = waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    if ~isempty(value) % check for error before executing switch
        switch value
            case 28 % left
                if current == 1
                    current = Moves;
                else
                    current = current - 1;
                end
            case 29 % right
                if current == Moves
                    current = 1;
                else
                    current = current + 1;
                end
            case 30 % up
                % do nothing for now
            case 31 % down
                go = 1;
            otherwise

        end
        
    end
    
end

decision = DATA.moves{pieceNum}{current};
DATAout = DATA;
