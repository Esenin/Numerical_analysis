% @arg c - coeff vector of polynom with zeros. First value - max degree,
% a_0 > 0
% @result bounds of roots
function [ downN, upperN, downP, upperP ] = rootBounds( c )
    upperP = upperPosBound(c);
% bottom bound of positive roots [
    c2 = fliplr(c);
    downP = 1 / upperPosBound(c2);
% ]
% bottom bound of negative roots [
    c3 = negativeTransposition(c);   
    
    downN = -upperPosBound(c3);

% ]
% upper bound of negative roots [
    c4 = negativeTransposition(c);
    c4 = fliplr(c4);
    upperN = -1 / upperPosBound(c4);
% ]

end

