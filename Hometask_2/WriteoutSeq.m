function [  ] = writeoutSeq( xks, yks, f, g )
%WriteoutSeq prints sequence of Newton algo results
%   @arg xks  sequence of Newton approximation
%   @arg yks -- same as xks
%   @arg f   @f(x,y) = @f(x,y,ai)
%   @arg g   -- same as f

    fprintf('i= \t|\tXn = \t\t|\tYn = \t\t|\tf(Xn,Yn) = \t|\tg(Xn,Yn) = \t\n');
    for i = 1 : length(xks)
        fr = f(xks(i), yks(i));
        gr = g(xks(i), yks(i));
        fprintf('%i\t| %i\t| %i\t| %i\t| %i\t\n', i, xks(i), yks(i), fr, gr);
    end;
    fprintf('---------------------------------------------------------------\n');


end

