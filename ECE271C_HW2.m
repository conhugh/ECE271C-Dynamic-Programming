  %Connor Hughes
%ECE 271C HW2

%% Problem 2: World Health Council Problem
%initialize network of states:
nodes = [0, 0, 0; 1, 1, 1; 2, 2, 2; 3, 3, 3; ...
         4, 4, 4; 5, 5, 5];
    
%initialize matrix of person lives saved (given in problem)
PLS = [0, 0, 0; 45, 20, 50; 70, 45, 70; 90, 75, 80; 105, 110, 100; 120, 150, 130];
%negate person lives saved matrix so we can "minimize cost" and find
%optimal solution
PLS = -PLS;

%generate cost-to-go matrix:
L = Inf(6, 3);
last_node = zeros(6, 3);
%initialize first stage costs
L(:, 1) = flip(PLS(:, 1));
%compute cost-to-go for each of the other states:
for k = 1:2
   for j = 1:6
       cost_opts = Inf(6, 1);
       for i = j:6
           diff = (nodes(i, k) - nodes(j, k + 1)) + 1;
           if diff > 0
               cost_opts(i) = L(i, k) + PLS(diff, k + 1);
           end
       end
       [L(j, k + 1), last_node(j, k + 1)] = min(cost_opts);
   end
end

%assemble optimal path (converting result to teams assigned at each stage) 
%and find min overall ctg:
[min_ctg, index] = min(L(:, 3));
diff2 = nodes(last_node(index, 3), 2) - nodes(index, 3);
diff1 = nodes(last_node(last_node(index, 3), 2), 1) - nodes(last_node(index, 3), 2);
diff0 = 5 - nodes(last_node(last_node(index, 3), 2), 1);
path = [diff0, diff1, diff2]
min_ctg
%% Problem 3: Government Space Project 
%initialize network of states:
nodes = [0, 0, 0; 1, 1, 1; 2, 2, 2];

%initialize matrix of failure probabilities:
PF = [0.4, 0.6, 0.8; 0.2, 0.4, 0.5; 0.15, 0.2, 0.3];
%normalize to get relative reduction in failure probabilities:
PFR = zeros(3);
PFR(:, 1) = PF(:, 1)./PF(1, 1);
PFR(:, 2) = PF(:, 2)./PF(1, 2);
PFR(:, 3) = PF(:, 3)./PF(1, 3);

%generate cost-to-go matrix:
L = Inf(3, 3);
last_node = zeros(3, 3);
%initialize first stage costs
L(:, 1) = log(flip(PFR(:, 1)));
%compute cost-to-go for each of the other states:
for k = 1:2
   for j = 1:3
       cost_opts = Inf(3, 1);
       for i = j:3
           diff = (nodes(i, k) - nodes(j, k + 1)) + 1;
           if diff > 0
               cost_opts(i) = L(i, k) + log(PFR(diff, k + 1));
           end
       end
       [L(j, k + 1), last_node(j, k + 1)] = min(cost_opts);
   end
end

%assemble optimal path (converting result to scientists assigned at each stage) 
%and find min overall ctg:
[min_ctg, index] = min(L(:, 3));
diff2 = nodes(last_node(index, 3), 2) - nodes(index, 3);
diff1 = nodes(last_node(last_node(index, 3), 2), 1) - nodes(last_node(index, 3), 2);
diff0 = 2 - nodes(last_node(last_node(index, 3), 2), 1);
path = [diff0, diff1, diff2]
min_ctg
min_prob = 0.192*exp(min_ctg)
%% Problem 4: Markov Chain
%initialize matrix of transition probabilities:
P = [0.7 0.3 0; 0.5 0.1 0.4; 0 0.6 0.4];
%initialize signal observtion probability matrix:
Z = [0.5 0.5 0; 1/3 1/3 1/3; 0 0.5 0.5];
%initialize sequence of signals:
sigs = [1, 2, 3, 2];

X_k = [0.5 0.5 0];
X_kplusone = (X_k*P).*(Z(:, sigs(1)))';
X_kplusone_check = (X_k*P).*((Z(:, sigs(1)))./sum(Z(:, sigs(1))))';
Z_check = ((Z(:, sigs(1)))./sum(Z(:, sigs(1))))';
Xkp_check = (X_k*P);
X_kplusone = X_kplusone./(sum(X_kplusone));  %normalize result

%Part i:
%initial beliefs:
q_0 = [0.5, 0.5, 0];
q_10 = q_0*(P^10);

%Part ii:
%generate random initial condition according to initial belief:
x = zeros(11, 1);
x(1) = randsample([1, 2, 3], 1, true, [0.5, 0.5, 0]);
results = zeros(1000, 1);
for s = 1:1000
    for k = 1:10
        x(k + 1) = randsample([1, 2, 3], 1, true, P(x(k), :));
    end
    results(s) = x(11);
end
[GC, GR] = groupcounts(results)
histogram(results)
xticks([1, 2, 3])
xlabel("State at Stage 10")
ylabel("Number of occurrences")
title("Markov Chain Simulation")
emp_freq = GC./sum(GC)

%Part iii:
%initialize beliefs matrix:
q = zeros(5, 3);
%populate initial beliefs:
q(1, :) = [0.5 0.5 0];
%iterate through each stage/transition:
for k = 1:4
    denom = 0;
    %for each node transitioned to at stage k + 1:
    for i = 1:3
       %for each node potentially transitioned from at stage k:
       for j = 1:3
           %add probability of the transition, given known signal, to the
           %probability of arriving at current node in stage k + 1
           q(k + 1, i) = q(k + 1, i) + Z(i, sigs(k))*P(j, i)*q(k, j);
           %add up the total probability of seing signal
           denom = denom + Z(i, sigs(k))*P(j, i)*q(k, j);
       end
    end
    %compute belief at stage k + 1
    q(k + 1, :) = q(k + 1, :)./denom;
end

%Part iv: Viterbi Algorithm
%generate matrix of costs-to-go:
costs = Inf(3, 5);
for n = 1:3
   costs(n, 1) = -log(q(1, n));
end
last_node = zeros(3, 5);
for k = 1:4
   for j = 1:3
       cost_opts = Inf(3, 1);
       for i = 1:3
           cost_opts(i) = costs(i, k) - log(P(i, j)*Z(j, sigs(k)));
       end
       [costs(j, k + 1), last_node(j, k + 1)] = min(cost_opts);
   end
end

%find min total cost-to-go and assemble path:
paths = zeros(3, 5);
paths(:, 5) = [1; 2; 3];
for k = 4:-1:1
   for n = 1:3
       paths(n, k) = last_node(paths(n, k + 1), k + 1);
   end
end
[min_ctg, path_ind] = min(costs(:, 5));
path = paths(path_ind, :);
%% Helper Functions
function combs = nmultichoosek(values, k)
    %// Return number of multisubsets or actual multisubsets.
    if numel(values)==1
        n = values;
        combs = nchoosek(n+k-1,k);
    else
        n = numel(values);
        combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
        combs = reshape(values(combs),[],k);
    end
end