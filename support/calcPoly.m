function [ y ] = calcPoly( coeff, x0 )
%    y = f(x0)
    y = 0;
    for i = 1 : length(coeff)
        degree = length(coeff) - i;
        y = y + coeff(i) *(x0^degree);
    end;
end

