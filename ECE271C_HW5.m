%Connor Hughes
%ECE 271C HW5

%% Problem 5.1: Tennis
%initialize parameters:
qf = 0.6;
qs = 0.4;
ps = 0.95;
pfOpts = linspace(0, 1, 101);
nStages = 50;

%initialize transition probability matrices:
Pf = zeros(8, 8);  
Ps = zeros(8, 8);

Ps(1, 1) = 1;
Ps(2, 1) = qs*ps;
Ps(2, 3) = 1 - ps;
Ps(2, 4) = (1 - qs)*ps;
Ps(3, 1) = qs*ps;
Ps(3, 4) = (1 - qs)*ps + (1 - ps);
Ps(4, 2) = qs*ps;
Ps(4, 5) = 1 - ps;
Ps(4, 6) = (1 - qs)*ps;
Ps(5, 2) = qs*ps;
Ps(5, 6) = (1 - qs)*ps + (1 - ps);
Ps(6, 4) = qs*ps;
Ps(6, 7) = 1 - ps;
Ps(6, 8) = (1 - qs)*ps;
Ps(7, 4) = qs*ps;
Ps(7, 8) = (1 - qs)*ps + (1 - ps);
Ps(8, 8) = 1;

RTGcomp = zeros(8, length(pfOpts));
OCCcomp = zeros(8, length(pfOpts));

for f = 1:length(pfOpts)
    pf = pfOpts(f);
    Pf(1, 1) = 1;
    Pf(2, 1) = qf*pf;
    Pf(2, 3) = 1 - pf;
    Pf(2, 4) = (1 - qf)*pf;
    Pf(3, 1) = qf*pf;
    Pf(3, 4) = (1 - qf)*pf + (1 - pf);
    Pf(4, 2) = qf*pf;
    Pf(4, 5) = 1 - pf;
    Pf(4, 6) = (1 - qf)*pf;
    Pf(5, 2) = qf*pf;
    Pf(5, 6) = (1 - qf)*pf + (1 - pf);
    Pf(6, 4) = qf*pf;
    Pf(6, 7) = 1 - pf;
    Pf(6, 8) = (1 - qf)*pf;
    Pf(7, 4) = qf*pf;
    Pf(7, 8) = (1 - qf)*pf + (1 - pf);
    Pf(8, 8) = 1;

    %initialize reward-to-go matrix:
    RTG = zeros(8, nStages + 1);
    RTG(:, nStages + 1) = [1; 0; 0; 0; 0; 0; 0; 0];

    %initialize control matrix:
    OCC = zeros(8, nStages);

    for k = nStages:-1:1
       rtgF = Pf*RTG(:, k + 1);
       rtgS = Ps*RTG(:, k + 1);
       [RTG(:, k), OCC(:, k)] = max([rtgF'; rtgS']);
    end
    RTGcomp(:, f) = RTG(:, 1);
    OCCcomp(:, f) = OCC(:, 1);
end

%plot results:
plot(pfOpts, RTGcomp(4, :));
xlabel('pf');
ylabel('P(Server wins)')
title('Games Beginning at 30-30')
ylim([0 1])

%% Problem 5.2: Football!
%initialize parameters for relevant probabilities:
params = struct;
params.lamR = 3;
params.lamP = 10;
params.p = 0.4;
params.q = 0.05;

%initialize number of states and stages:
nStates = 1202;  %1-30 yards to goal * 1-10 yards to first * 1-4 downs + 2 gameover
nStages = 12;

%initialize network of states:
states = cell(nStates, nStages + 1);
for k = 1:nStages + 1
   for i = 1:(nStates - 2)
       %populate states with the following format:
       % [yards to goal, down, yards to first]
       states{i, k} = zeros(1, 3);
       states{i, k}(1) = floor((i - 1)/40) + 1;  %initialize yards to goal
       states{i, k}(2) = floor(mod((i - 1), 40)/10) + 1;
       states{i, k}(3) = mod(i - 1, 10) + 1;
   end
   states{nStates - 1, k} = "TD";  %initialize touchdown state
   states{nStates, k} = "EOD";  %initialize end of drive (turnover) state
end

%initialize matrix of "next nodes":
NN = NaN(nStates, nStages); 

%initialize matrix of optimal control choices:
OCC = NaN(nStates, nStages);

%initialize reward-to-go matrix:
RTG = NaN(nStates, nStages + 1);
RTG(:, nStages + 1) = zeros(nStates, 1);
RTG(nStates - 1, nStages + 1) = 1;

%populate the reward-to-go matrix:
for k = nStages:-1:1
    for i = 1:nStates    %for each state at each stage
        pTest = zeros(nStates, 1);
        rTest = zeros(nStates, 1);
        rtgR = 0;
        rtgP = 0;
        %compute expected reward-to-go for running and for passing
        for j = 1:nStates
            rtgR = rtgR + transProb(states{i, k}, 0, states{j, k + 1}, params)*RTG(j, k + 1);  
            rtgP = rtgP + transProb(states{i, k}, 1, states{j, k + 1}, params)*RTG(j, k + 1);
        end
        %set reward-to-go equal to the max of the two options, set optimal 
        %control choice equal to the choice which maximized reward-to-go:
        [RTG(i, k), OCC(i, k)] = max([rtgR, rtgP]);   
    end
end

% Plot Results
close all;
%first down:
ind = 1:40:1200;
vC = zeros(1, 10*length(ind));
for i = 1:length(ind)
    j = 10*(i - 1) + 1;
    vC(j:(j + 9)) = OCC(ind(i):(ind(i) + 9), 1);
end
mC = reshape(vC, [10, 30]);
figure
heatmap(mC)
title('First Down')
xlabel('Yards to End Zone')
ylabel('Yards to First Down')

%second down:
ind = 11:40:1200;
vC = zeros(1, 10*length(ind));
for i = 1:length(ind)
    j = 10*(i - 1) + 1;
    vC(j:(j + 9)) = OCC(ind(i):(ind(i) + 9), 1);
end
mC = reshape(vC, [10, 30]);
figure
heatmap(mC)
title('Second Down')
xlabel('Yards to End Zone')
ylabel('Yards to First Down')

%third down:
ind = 21:40:1200;
vC = zeros(1, 10*length(ind));
for i = 1:length(ind)
    j = 10*(i - 1) + 1;
    vC(j:(j + 9)) = OCC(ind(i):(ind(i) + 9), 1);
end
mC = reshape(vC, [10, 30]);
figure
heatmap(mC)
title('Third Down')
xlabel('Yards to End Zone')
ylabel('Yards to First Down')

%fourth down:
ind = 31:40:1200;
vC = zeros(1, 10*length(ind));
for i = 1:length(ind)
    j = 10*(i - 1) + 1;
    vC(j:(j + 9)) = OCC(ind(i):(ind(i) + 9), 1);
end
mC = reshape(vC, [10, 30]);
figure
heatmap(mC)
title('Fourth Down')
xlabel('Yards to End Zone')
ylabel('Yards to First Down')

%% Problem 5.3: Computer Manufacturer
%initialize parameters:
p11A = 0.8;
p12A = 0.2;
p11N = 0.5;
p12N = 0.5;
p21R = 0.7;
p22R = 0.3;
p21N = 0.4;
p22N = 0.6;
alpha = 0.99;
nStages = 1000;

%initialize reward-to-go and control matrices:
RTG = zeros(2, nStages + 1);
OCC = zeros(2, nStages);

cUpr = NaN;
cLwr = NaN;
uprBnd = NaN(2, nStages);
lwrBnd = NaN(2, nStages);

%populate reward-to-go matrix:
for k = nStages:-1:1
   for i = 1:2
       if i == 1
           rtgA = 4 + alpha*(p11A*RTG(1, k + 1) + p12A*RTG(2, k + 1));
           rtgN = 6 + alpha*(p11N*RTG(1, k + 1) + p12N*RTG(2, k + 1));
           rtgOpts = [rtgA, rtgN];
       else
          rtgR = -5 + alpha*(p21R*RTG(1, k + 1) + p22R*RTG(2, k + 1));
          rtgN = -3 + alpha*(p21N*RTG(1, k + 1) + p22N*RTG(2, k + 1));
          rtgOpts = [rtgR, rtgN];
       end
       [RTG(i, k), ind] = max(rtgOpts);
       diff = (RTG(:, k) - RTG(:, k + 1));
       cUpr = max(diff); 
       cLwr = min(diff); 
       uprBnd(:, k + 1) = (alpha/(1 - alpha))*cUpr + RTG(:, k);
       lwrBnd(:, k + 1) = (alpha/(1 - alpha))*cLwr + RTG(:, k);
       OCC(i, k) = (i - 1)*2 + ind;
   end
end

%plot results
figure
stages = linspace(1, nStages, nStages);
plot(stages, RTG(1, 1:(end - 1)), stages, RTG(2, 1:(end - 1)), stages, uprBnd(1, 1:(end - 1)), stages, uprBnd(2, 1:(end - 1)), stages, lwrBnd(1, 1:(end - 1)), stages, lwrBnd(2, 1:(end - 1)), 'LineWidth', 1.2)
legend('RTG from State 1', 'RTG from State 2', 'Upper Bound on RTG from State 1', 'Upper Bound on RTG from State 2', 'Lower Bound on RTG from State 1', 'Lower Bound on RTG from State 2')
title('Reward-to-go Convergence')
xlabel('Stage Number')
ylabel('Reward-to-go')
ax = gca
ax.FontSize = 20

figure
stages = linspace(1, 10, 10);
plot(stages, uprBnd(1, (end - 10):(end - 1)), stages, uprBnd(2, (end - 10):(end - 1)), stages, lwrBnd(1, (end - 10):(end - 1)), stages, lwrBnd(2, (end - 10):(end - 1)), 'LineWidth', 1.2)
legend('Upper Bound on RTG from State 1', 'Upper Bound on RTG from State 2', 'Lower Bound on RTG from State 1', 'Lower Bound on RTG from State 2')
title('Reward-to-go Convergence (Error Bounds)')
xlabel('Stage Number')
ylabel('Reward-to-go')
ax = gca
ax.FontSize = 20

%optimal policy comparison for alpha = 0.9:
a = 0.9;

c1 = 4;
c2 = -5;
p11 = 0.8;
p22 = 0.3;
J1 = 1/(1 - a*(p11 + p22) - a^2*(1 - p11 - p22))*(c1 - (c1*p22 - c2 + c2*p11)*a)
J2 = 1/(1 - a*(p11 + p22) - a^2*(1 - p11 - p22))*(c2 + (c1 - c1*p22 - c2*p11)*a)

c1 = 6;
c2 = -5;
p11 = 0.5;
p22 = 0.3;
J1 = 1/(1 - a*(p11 + p22) - a^2*(1 - p11 - p22))*(c1 - (c1*p22 - c2 + c2*p11)*a)
J2 = 1/(1 - a*(p11 + p22) - a^2*(1 - p11 - p22))*(c2 + (c1 - c1*p22 - c2*p11)*a)

c1 = 4;
c2 = -3;
p11 = 0.8;
p22 = 0.6;
J1 = 1/(1 - a*(p11 + p22) - a^2*(1 - p11 - p22))*(c1 - (c1*p22 - c2 + c2*p11)*a)
J2 = 1/(1 - a*(p11 + p22) - a^2*(1 - p11 - p22))*(c2 + (c1 - c1*p22 - c2*p11)*a)

c1 = 6;
c2 = -3;
p11 = 0.5;
p22 = 0.6;
J1 = 1/(1 - a*(p11 + p22) - a^2*(1 - p11 - p22))*(c1 - (c1*p22 - c2 + c2*p11)*a)
J2 = 1/(1 - a*(p11 + p22) - a^2*(1 - p11 - p22))*(c2 + (c1 - c1*p22 - c2*p11)*a)

%% Problem 5.5: Umbrella Problem
%initialize parameters:
syms p alpha W V 
assume(p >= 0);
assume(p <= 1);
assume(alpha > 0);
assume(alpha < 1);
assume(W >= 0);
assume(V >= 0);

%define transition probabilities:
Pt = [p, 0, 1-p, 0; p, 0, 1-p, 0; p, 0, 1-p, 0; p, 0, 1-p, 0];
Pnt = [p, 0, 1-p, 0; p, 0, 1-p, 0; 0, p, 0, 1-p; p, 0, 1-p, 0];

%define stage costs:
Gt = [0; W; V; 0];
Gnt = [0; W; 0; 0];

%Mu-Specific Bellman Equations:
Jt = inv(eye(4) - alpha*Pt)*Gt;
Jnt = inv(eye(4) - alpha*Pnt)*Gnt;

cond = Jt <= Jnt;

simplify(cond(3))
%% Helper Functions:
%Football problem: returns probability of transition from state
%xk to xkp, given control uk (1 = pass, 0 = run)
function p = transProb(xk, uk, xkp, params)
    if xk(1) == 'TD'       %if current state is touchdown
        if xkp(1) == 'TD'  %if next state is touchdown
           p = 1;
        else
           p = 0;
        end
    elseif xk(1) == 'EOD'       %if current state is end-of-drive
        if xkp(1) == 'EOD'  %if next state is end-of-drive
           p = 1;
        else
           p = 0;
        end
    elseif xkp(1) == 'TD'   %if next state is touchdown
        if uk == 0      %if the control choice is to run
            p = 1 - poisscdf(xk(1) - 1, params.lamR);
        end
        if uk == 1      %if the control choice is to pass
            p = (1 - params.p - params.q)*(1 - poisscdf(xk(1) - 1, params.lamP));
        end
    elseif xkp(1) == 'EOD'   %if next state is end of drive
        if uk == 0       %if the control choice is to run
            if xk(2) ~= 4   %if it's not 4th down
               p = 0;
            else
               p = poisscdf(min([xk(1), xk(3)]) - 1, params.lamR);  %return prob(run was short of 1st down) 
            end
        end
        if uk == 1       %if the control choice is to pass
            if xk(2) ~= 4   %if it's not 4th down
               p = params.q;   %return prob(pass intercepted)
            else
               %return prob(pass intercepted OR incomplete OR short of 1st down/goal line): 
               p = params.q + params.p + (1 - params.q - params.p)*poisscdf(min([xk(1), xk(3)]) - 1, params.lamP);   
            end
        end
    else
        %next state is neither TD nor EOD
        yg = -(xkp(1) - xk(1));  %compute yards gained
        %if either down number or number of yards to 1st down
        %are incorrect at the next stage:
        if (xkp(2) ~= nextDown(xk, yg)) || (xkp(3) ~= nextYTF(xk, yg))
            p = 0;
        else
            if uk == 0        %if the control choice is to run
                p = poisspdf(yg, params.lamR);
            end
            if uk == 1        %if the control choice is to pass
                if yg == 0    %if there are zero yards gained
                    p = params.p + (1 - params.p - params.q)*poisspdf(yg, params.lamP);
                else
                    p = (1 - params.p - params.q)*poisspdf(yg, params.lamP);
                end
            end
        end
    end
    
end

%Football problem: returns the appropriate next down number, given 
%current state xk and number of yards gained yg (assuming neither a TD 
%or TO occurs)
function d = nextDown(xk, yg)
    if yg >= xk(3)      %if enough yards gained for 1st down
        d = 1;
    else
        d = xk(2) + 1; 
    end
end

%Football problem: returns the appropriate next number of yards-to-1st-down, 
%given current state xk and number of yards gained yg (assuming neither a TD 
%or TO occurs)
function ytf = nextYTF(xk, yg)
    if yg >= xk(3)       %if enough yards gained for 1st down
       ytf = 10;
    else
       ytf = xk(3) - yg; 
    end
end