%Connor Hughes
%ECE 271C HW1

%% Tests
%Problem 5:
L = zeros(3, 3, 2);
L(:, :, 1) = [10 20 2; 1 2 3; 2 1 0];
L(:, :, 2) = [0 1 2; 1 0 2; 2 2 1];
route = pathplan(L)

%Problem 6
X = [7 9 1 4 13 2];
comps_a = matrix_multiply_a(X)
comps_b = matrix_multiply_b(X)

%% Problem 5:
function route = pathplan(L)
    %initialize cost-to-go vector for stage k:
    ctg_k = zeros(length(L(:, 1, 1)), 1);
    %initialize cost-to-go vector for stage k+1:
    ctg_kplusone = zeros(length(L(:, 1, 1)), 1);
    %initialize vector for storing cost options from a given node:
    total_costs = zeros(length(L(:, 1, 1)), 1);
    %initialize vectors for storing optimal paths from each node at stages k and k+1
    paths = strings(length(L(:, 1, 1)), 1);
    paths_kplusone = strings(length(L(:, 1, 1)), 1);
    paths_kplusone = paths_kplusone + "t";
    %iterate backward through each stage:
    for k = (length(L(1, 1, :))):-1:1
        %iterate through each node at current stage:
        for i = 1:length(L(:, 1, 1))
            %iterate through each target node for the current node:
            for j = 1:length(L(1, :, 1))
                %evaluate and store costs for tail problem:
                total_costs(j) = L(i, j, k) + ctg_kplusone(j);
            end
            %find minimum cost for tail problem:
            [min_cost, index] = min(total_costs);
            %store min cost as cost-to go from this node, and store
            %corresponding path
            ctg_k(i) = min_cost;
            paths(i) = i + paths_kplusone(index);
        end
        ctg_kplusone = ctg_k;
        paths_kplusone = paths;
    end
    %evaluate min cost-to-go from starting node, and return path:
    [min_cost, index] = min(ctg_k);
    route = "s" + paths(index);
end

%% Problem 6:
%Part A:
function computations_a = matrix_multiply_a(X)
    %enumerate the multiplications:
    mults = [1:1:(length(X) - 2)];
    %generate the network of nodes/states:
    num_stgs = length(mults) + 1;
    nodes = cell(1, num_stgs);
    for k = 1:num_stgs
        temp = nchoosek(mults, num_stgs - k);
        nodes{k} = temp;
    end

    %generate cell structure with matrices of stage costs, initialized to zeros:
    L = cell(1, num_stgs);
    for k = 1:(num_stgs - 2) 
        L{k} = zeros(length(nodes{k}(:, 1)), length(nodes{k+1}(:, 1)));
    end
    L{num_stgs - 1} = zeros(length(nodes{num_stgs - 1}), 1);
    
    %for each stage except the termination stage
    for k = 1:(num_stgs - 1)
       %for each node in the stage
       for i = 1:length(L{k}(:, 1)) 
           %for each target node in the next stage
            for j = 1:length(L{k}(1, :)) 
                %for each multiplication remaining in the target node at the next stage
                for p = 1:length(nodes{k + 1}(1, :))
                    %if it's not in the current node, set cost to infinite
                    if ~ismember(nodes{k+1}(j, p), nodes{k}(i, :))
                        L{k}(i, j) = Inf;
                    end
                end
                %if the cost to go from the current node to the target node is
                %not infinite
                if L{k}(i, j) ~= Inf
                    %for each multiplication in the current node
                    for q = 1:length(nodes{k}(1, :))
                        %if it's not listed in the target node
                        if ~ismember(nodes{k}(i, q), nodes{k+1}(j, :))
                            %conclude that it is the multiplication which was
                            %completed, determine number of computations, and
                            %store this value as the cost to get from the
                            %current node to target node
                            u = nodes{k}(i, q);
                            state = nodes{k}(i, :);
                            L{k}(i, j) = X(max([0, state(state<u)]) + 1)*X(u+1)*X(min([num_stgs, state(state>u)]) + 1);                                
                        end
                    end
                end
            end
       end    
    end
    %Now, with cost matrix populated, use reverse DP algorithm to solve:
    ctg_kplusone = zeros(1);
    for k = (length(L) - 1):-1:1
        ctg_k = Inf(length(L{k}(:, 1)), 1);
        for i = 1:length(L{k}(:, 1))
            total_costs = zeros(length(L{k}(1, :)), 1);
            for j = 1:length(L{k}(1, :))
                total_costs(j) = L{k}(i, j) + ctg_kplusone(j);
            end
            min_cost = min(total_costs);
            ctg_k(i) = min_cost;
        end
        ctg_kplusone = ctg_k;
    end
    computations_a = min(ctg_k);
end
%% Problem 6
%Part B:
function computations_b = matrix_multiply_b(X)
    %enumerate the multiplications:
    mults = [1:1:(length(X) - 2)];
    %generate the network of nodes/states:
    num_stgs = length(mults) + 1;
    nodes = cell(1, num_stgs);
    for k = 1:num_stgs
        temp = nchoosek(mults, num_stgs - k);
        nodes{k} = temp;
    end

    %generate cell structure with matrices of stage costs, initialized to zeros:
    L = cell(1, num_stgs);
    for k = 1:(num_stgs - 2) 
        L{k} = zeros(length(nodes{k}(:, 1)), length(nodes{k+1}(:, 1)));
    end
    L{num_stgs - 1} = zeros(length(nodes{num_stgs - 1}), 1);
    
    %for each stage except the termination stage
    for k = 1:(num_stgs - 1)
       %for each node in the stage
       for i = 1:length(L{k}(:, 1)) 
           %for each target node in the next stage
            for j = 1:length(L{k}(1, :)) 
                %for each multiplication remaining in the target node at the next stage
                for p = 1:length(nodes{k + 1}(1, :))
                    %if it's not in the current node, set cost to infinite
                    if ~ismember(nodes{k+1}(j, p), nodes{k}(i, :))
                        L{k}(i, j) = Inf;
                    end
                end
                %if the cost to go from the current node to the target node is
                %not infinite
                if L{k}(i, j) ~= Inf
                    %for each multiplication in the current node
                    for q = 1:length(nodes{k}(1, :))
                        %if it's not listed in the target node
                        if ~ismember(nodes{k}(i, q), nodes{k+1}(j, :))
                            %conclude that it is the multiplication which was
                            %completed, determine number of computations, and
                            %store this value as the cost to get from the
                            %current node to target node
                            u = nodes{k}(i, q);
                            state = nodes{k}(i, :);
                            %check if this multiplication was valid (if it 
                            %was adjacent to those already completed)
                            if (k > 1 && ((ismember(u + 1, state) && ismember(u - 1, state)) ...
                                   || (u == 1 && ismember(u + 1, state)) || (u == num_stgs - 1 && ismember(u - 1, state))))
                                %if not, set cost to inf
                                L{k}(i, j) = Inf;
                            else
                                %evaluate computational cost
                                L{k}(i, j) = X(max([0, state(state<u)]) + 1)*X(u+1)*X(min([num_stgs, state(state>u)]) + 1);                                
                            end
                        end
                    end
                end
            end
       end    
    end
    %Now, with cost matrix populated, use reverse DP algorithm to solve:
    ctg_kplusone = zeros(1);
    for k = (length(L) - 1):-1:1
        ctg_k = Inf(length(L{k}(:, 1)), 1);
        for i = 1:length(L{k}(:, 1))
            total_costs = zeros(length(L{k}(1, :)), 1);
            for j = 1:length(L{k}(1, :))
                total_costs(j) = L{k}(i, j) + ctg_kplusone(j);
            end
            min_cost = min(total_costs);
            ctg_k(i) = min_cost;
        end
        ctg_kplusone = ctg_k;
    end
    computations_b = min(ctg_k);
end
