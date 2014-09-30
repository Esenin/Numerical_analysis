clear 
clc

%imports
addpath('../support/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    9x4 -9x2 -36x + 1
polyExample     = @(x) 9*x^4 - 9*x^2 - 36*x + 1;
polyExampleAprx = [9, 0, -9, -36, 1];
polyExampleDer  =  @(x) 36*x^3 -18*x -36;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    x3 - 4x + 2 = 0
poly11     = @(x) x^3 - 4*x + 2;
poly11Aprx = [1 0 -4 2];
poly11Der  = @(x) 3*x^2 - 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    x3 + sin x - 12x + 1
poly12     = @(x) x^3 + sin(x) - 12*x + 1;
poly12Aprx = [(-121*10^(-15)) 0 (355*10^(-12)) 0 (-1.3*10^(-12)) 0 (6.2*10^(-9)) 0 (-3.9*10^(-7)) 0 (1/362880) 0 (-1/5040) 0 (1/120) 0 (5/6) 0 -11 1];
poly12Der  = @(x) 3*x^2 + cos(x) - 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    3x6 + 172x5 + 125x4 -600x3 +125x2 + 172x + 3 
poly6     = @(x) 3*x^6 + 172*x^5 + 125*x^4 -600*x^3 +125*x^2 + 172*x + 3;
poly6Aprx = [3 172 125 -600 125 172 3];
poly6Der  = @(x) 18*x^5 + 860*x^4 + 500*x^3 -1800*x^2 +250*x + 172;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long;
tic

%egRoots = NewtonRoots(polyExample, polyExampleAprx, polyExampleDer)

%eqOneRoots = NewtonRoots(poly11, poly11Aprx, poly11Der)

poly12
eqTwoRoots = NewtonRoots(poly12, poly12Aprx, poly12Der)

%eq6Roots = NewtonRoots(poly6, poly6Aprx, poly6Der)

toc


