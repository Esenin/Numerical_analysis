function [ nc ] = negativeTransposition( coeff )
%   transposition: x = -x;
    nc = coeff;
    for i = 1 : length(nc)
        curDeg = length(nc) - i;
        if (mod(curDeg, 2) == 1)
            nc(i) = -nc(i);
        end;
    end;

end

