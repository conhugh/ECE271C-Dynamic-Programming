function [decision,DATAout] = myPlay2(~,pieceNum,DATA)
% I think this is deterministic?

global STAGENUMBER;
stage = STAGENUMBER;
load 'control.mat' control;
decision = DATA.moves{pieceNum}{control(stage)};
DATAout = DATA;
STAGENUMBER = STAGENUMBER + 1;
