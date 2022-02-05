function [decision,DATAout] = TDLearning(board,pieceNum,DATA)
% Testing out various policies
% 
% Places the piece to minimize the maximum height
%

% global WEIGHTS;
% R = WEIGHTS
global TIME;
T = TIME;
global ZVEC;
Z = ZVEC;
global WEIGHTUPDATE;
RUPD = WEIGHTUPDATE;
R = RUPD

% I wonder if some sort of diagonalization/centering/scaling would be good?

cur_phi = extract_features(board);
%cur_est_ctg = R*cur_phi; % using updated weights to estimate

next_phi = NaN(2,length(DATA.moves{pieceNum}));
scores = NaN(1,length(DATA.moves{pieceNum}));
%u_est_ctg = NaN(1,length(DATA.moves{pieceNum}));

for u = 1:length(DATA.moves{pieceNum})
    
    move = DATA.moves{pieceNum}{u};
    [theNextBoard,nextScore] = nextBoard(board,move);
    
    next_phi(:,u) = extract_features(theNextBoard);
    scores(u) = nextScore;
    
end

%CONNOR ADDED THIS TO CHECK HOW SCORES ARE COUNTED:
% if(max(scores) ~= 0)
%    qrpfh = 0; %useless garbage variable, just here to set a breakpoint
% end

% Estimate cost-to-go for stable weights
alpha = 1;
u_est_ctg = alpha*R*next_phi - scores

% Find best control and choose via softmax
[~,u_best] = min(u_est_ctg)
u_aug = zeros(1,length(u_est_ctg));
u_aug(u_best) = 1;
u_smax = smax(u_aug,0.2) % change T parameter?
u_act = randsample(1:length(u_smax),1,true,u_smax)

%orig_Tdiff = alpha*R*next_phi(:,u_best) - scores(u_best) - cur_imp_ctg

% Improve estimate of weights
cur_imp_ctg = RUPD*cur_phi;
imp_est_ctg = alpha*RUPD*next_phi(:,u_best) - scores(u_best);
Tdiff = imp_est_ctg - cur_imp_ctg
gamma = 1/T;
lambda = 0;
Zt = lambda*Z + cur_phi;

next_weights = RUPD - gamma*Tdiff*Zt';
next_imp_ctg = next_weights*cur_phi;
next_Tdiff = alpha*next_weights*next_phi(:,u_best) - scores(u_best) - next_imp_ctg

WEIGHTUPDATE = RUPD - gamma*Tdiff*Zt';
ZVEC = Zt;
TIME = TIME + 1;

%decision = DATA.moves{pieceNum}{u_act};
decision = DATA.moves{pieceNum}{u_act};
DATAout = DATA;