function output = smax(input,t)
    %input = input/sum(input);
    output = exp(input/t)/sum(exp(input/t));
end

