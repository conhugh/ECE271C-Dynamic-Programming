function [decision,DATAout] = myPlay3(board,pieceNum,DATA)
% I think this is stochastic -- not sure if minimax or expected

global STAGENUMBER;
stage = STAGENUMBER;
load 'best_control.mat' best_control;

boardNum = b2d(flip(reshape(board,1,numel(board))))+1;
stateNum = (pieceNum-1)*512+boardNum;

control = best_control(stateNum,stage);

decision = DATA.moves{pieceNum}{control};
DATAout = DATA;
STAGENUMBER = STAGENUMBER + 1;

