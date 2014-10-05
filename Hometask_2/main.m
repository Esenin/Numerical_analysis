clear 
clc

%imports
addpath('../support/');
format long;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  { tg(ax + y) - axy = 0.3
%  { x^2 + y^2 = k
%    a = -0.6 : -0.1 : -1.1
%    k = 1.3 : 0.1 : 1.8
%a = -0.6 : -0.1 : -1.1;
%k = 1.3 : 0.1 : 1.8;
a = -1.1;
k = 1.8;

f = @(x, y, a) tan(a*x + y) - a*x*y - 0.3;
g = @(x, y, k) x^2 + y^2 - k;

fdx = @(x, y, a) 2 * a / (cos(2 * a * x + 2 * y) + 1)  - a * y ;
fdy = @(x, y, a) 2 / (cos(2 * a * x + 2 * y) + 1) - a * x;
gdx = @(x, ~, ~) 2 * x; 
gdy = @(~, y, ~) 2 * y; 


tic
NewtonRoots(f, g, fdx, fdy, gdx, gdy, a, k);
toc

figure
hold on
grid on
ezplot('tan(-0.6*x + y) + 0.6*x*y - 0.3',[-2 2 -2 2]);
ezplot('x^2 + y^2 - 1.8',[-2 2 -2 2]);

legend('tan(-0.6*x + y) + 0.6*x*y - 0.3',  ...
       'x^2 + y^2 - 1.8');



