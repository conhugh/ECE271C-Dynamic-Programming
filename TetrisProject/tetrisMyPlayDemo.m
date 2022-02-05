% Tetris demo 1
function [Iscore,nPieces] = tetrisMyPlayDemo(sData)

%     myFun = 'myPlay';
%    myFun = 'IPlay';
%    myFun = 'myPlay2';
%    myFun = 'myPlay3';
    myFun = 'TDLearning';
    global STAGENUMBER;
    STAGENUMBER = 1;
    
    global ZVEC;
    global WEIGHTS;
    global WEIGHTUPDATE;
    global TRACKWEIGHTS;
    global TIME;

    
    DATA = sData;
    
    startPiece = DATA.startPiece;
    GameSize = DATA.GameSize;
    N = DATA.nMaxPieces;
    RowCap = DATA.RowCap;
    Pcolor = DATA.Pcolor;
    Pieces = DATA.Pieces;
    nEpisodes = DATA.nEpisodes;
    S_Plot = DATA.S_Plot;
    TimeDelay = DATA.TimeDelay;
    S_Sounds = DATA.S_Sounds;
    Deterministic = DATA.Deterministic;
    
    for kc=1:nEpisodes
        
        % Tracks past basis vectors, reinitializes for each new game
        ZVEC = [0;0];
        WEIGHTS = WEIGHTUPDATE;
        TRACKWEIGHTS(:,TIME-100+1) = WEIGHTS';
        TIME = TIME + 1;
    
        CurPnum = startPiece;
        board=zeros(GameSize);

        if S_Plot==1
            figure(1)
        end
        GameOver=0;
        Score=0;
        stoploop=0;
    
        for numPlays = 1:N+1,
            if GameOver ~= 1,
            
                if Deterministic == 0
                    CurPnum = randi(length(Pieces));
                else
                    CurPnum = Deterministic;
                    if Deterministic == 3
                        Deterministic = 1;
                    else
                        Deterministic = Deterministic + 1;
                    end
                end

                % the game board that is under the red line
                boardLower = board(GameSize-RowCap+1:end,:); 
                boardLower = boardLower>0;
            
%                 theString = strcat('[decision,DATA] =',myFun,...
%                     '(boardLower,CurPnum,DATA);');
%                 eval(theString);
                [decision,DATA] = eval(strcat(myFun,'(boardLower,CurPnum,DATA);'));

                cols = find(sum(decision,1)>0);
                rows = find(sum(decision,2)>0);
                CurPlace = cols(1);
                CurPiece = decision(rows,cols);
                cp = Pieces{CurPnum};

                CurPiece = CurPiece*Pcolor(CurPnum);

                [CPr,CPc]=size(CurPiece);

                if CurPlace+CPc-1>GameSize(2)
                    warning('WARNING: Your piece is off the gameboard!')
                end

                overlap=0;
                GameBoard1=board;
                count=0;
                while overlap==0
                    count=count+1;

                    if count+CPr-2==GameSize(1)
                        overlap=1;
                        board=GameBoard1;
                    elseif max(max((board(count:count+CPr-1,CurPlace:CurPlace+CPc-1)>0)+(CurPiece>0)))>1
                        overlap=1;
                        board=GameBoard1;
                    else
                        GameBoard1=board;
                        GameBoard1(count:count+CPr-1,CurPlace:CurPlace+CPc-1)=board(count:count+CPr-1,CurPlace:CurPlace+CPc-1)+CurPiece;
                        %break;
                    end

                    if S_Plot==1

                        [r,c] = size(GameBoard1);                           % Get the matrix size
                        imagesc((1:c)+0.5,(1:r)+0.5,GameBoard1,[0,max(Pcolor)]);            % Plot the image
                        axis equal                                   % Make axes grid sizes equal
                        set(gca,'XTick',1:c+1,'YTick',1:r+1,...  % Change some axes properties
                        'XLim',[1 c+1],'YLim',[1 r+1],...
                        'XTickLabel',0:c,'YTickLabel',0:r,...
                        'GridLineStyle','-','XGrid','on','YGrid','on');
                        hold on
                        plot([0,r],[r-RowCap+1,r-RowCap+1],'r', 'LineWidth',3)
                        title(['Score=',num2str(Score), ', N=', num2str(numPlays)])
                        hold off
                        pause(TimeDelay);
                    end
                end
                
                nPieces(kc)=numPlays;
                S=sum(board>0,2)==GameSize(2);

                board(S==1,:)=[];
                Score=Score+sum(S);
                board=[zeros(sum(S),GameSize(2));board];
                if sum(S)>0 && S_Sounds==1
                    beep
                end

                if sum(sum(board(1:GameSize(1)-RowCap,:)))>0 || numPlays == N
                    GameOver=1;
                    if  S_Plot==1
                        title(['ROUND ', num2str(kc),...
                            ': Score=',num2str(Score),', Pieces=', num2str(numPlays)])
                        disp(['ROUND ', num2str(kc),...
                            ': Score=',num2str(Score),', Pieces=', num2str(numPlays)])
                    end
                    if S_Sounds==1, load gong.mat; sound(y,1*Fs);  end
                end
            end
            if GameOver==1,stoploop=1;  end
        end

        Iscore(kc)=Score; % Store the scores
    %fprintf('Score = %i\n',Score)
    
    end

end