function i = randiP(p)
%
% i = randiP(p)
%
% Selects a random integer according to probability vector p.
%
n = length(p);
p = cumsum(p);

t = rand;
for j = 1:n,
    if t <= p(j),
        i = j;
        return
    end
end

