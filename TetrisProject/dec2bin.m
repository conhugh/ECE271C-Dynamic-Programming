function [bin] = dec2bin(dec,nCol)
%DEC2BIN Summary of this function goes here
%   Detailed explanation goes here
    tmp = dec;
    bin = zeros(1,nCol);
    for i = 1:nCol
        bin(i) = mod(tmp,2);
        tmp = (tmp-bin(i))/2;
    end
end

