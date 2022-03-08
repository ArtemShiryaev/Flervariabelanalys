%-----------------------------------------------------------%   
%          Multi Variable Calculus (7.5 credits)            %
%          Olow Sande    Umeå University    2022            %
%-----------------------------------------------------------%   



clear % This line clears all saved variables and figures.


%................... TASKS .....................%


%...............................................%
%................. Functions ...................%
%...............................................%

% Scroll down to find an outline for the functions
% MyNPderiv, MyNewton and MyCritical. You should start by
% writing these functions.



% test1 = MyNPderiv( @(x,y) (x-1).^2 + (y-1).^2, 0, 0, 1)
% [test2x, test2y] = MyNewton( @(x,y) y+x , @(x,y) y-2 , 0, 0)
% [test3x, test3y] = MyCritical( @(x,y) (x-1).^2 + (y-1).^2, 0, 0)




%...............................................%
%.................. Figure 1 ...................%
%...............................................%

% Plot the graph of the topography function using surf.



% xline=linspace(0,10,100); % A vector with x-values from 0 to 10
% yline=linspace(0,10,100); % A vector with y-values from 0 to 10
% [x, y]=meshgrid(xline,yline); % Matrices of x and y values
% z = topography(x, y); % Calculates a matrix of z-values
% figure(1)
% surf(x, y, z, 'EdgeColor', 'none', 'FaceColor','interp')
% colormap summer




%...............................................%
%.................. Figure 2 ...................%
%...............................................%

% Plot the level curves of the topohraphy function using
% contour.




% figure(2)
% contour(x, y, z, 'g') % Level curves f(x,y)=C for topography




%...............................................%
%.................. Figure 3 ...................%
%...............................................%

% Plot 0-level curves of the partial derivatives (of
% topography) using contour.



% zx= MyNPderiv(@(x,y) topography(x,y), x, y, 1);
% zy= MyNPderiv(@(x,y) topography(x,y), x, y, 2);
% figure(3);
% contour(x, y, zx, [0,0], 'r')% Curves where df/dx=0
% hold on
% contour(x, y, zy, [0,0], 'p')% Curves where df/dy=0
% hold off




%...............................................%
%.................. Figure 4 ...................%
%...............................................%

% Plot the level curves of topography and the 0-level curves
% of the partial derivatives, in the same plot.



% figure(4);
% contour(x, y, z, 'g')% Level curves f(x,y)=C for topography
% hold on
% contour(x, y, zx, [0,0], 'r')% Curves where df/dx=0
% contour(x, y, zy, [0,0], 'p')% Curves where df/dx=0
% hold off




%...............................................%
%......... Minimum and maximum points ..........%
%...............................................%

% Use MyCritical to locate all local minimum and maximum
% points of the topography function. Remember that not every
% critical point is a local minimum or maximum point!



% [MaxX(1),MaxY(1)] = MyCritical(... , ... , ...)
% [MaxX(2),MaxY(2)] = MyCritical(... , ... , ...)
% ...
% [MinX(1),MinY(1)] = MyCritical(... , ... , ...)
% [MinX(2),MinY(2)] = MyCritical(... , ... , ...)
% ...





%...............................................%
%.................. Figure 5 ...................%
%...............................................%

% Plot the level curves, local maximum points and local
% minimum points of topography, in the same plot.



% figure(5);
% contour(x, y, z, 'g')% Level curves f(x,y)=C for topography
% hold on
% plot(MaxX, MaxY, 'r*')% Maximum points
% plot(MinX, MinY, 'b*')% Minimum points
% hold off




%...............................................%
%................ Figure 6 & 7 .................%
%...............................................%

% Implement the method Gradient Descent and plot a brook 
% running down the mountain.



% figure(6);
% contour(x, y, z, 'g')% Level curves f(x,y)=C for topography
% hold on
% plot(xgdes, ygdes, 'b')% A mountain brook!
% hold off

% figure(7);
% surf(x, y, z, 'EdgeColor', 'none', 'FaceColor','interp')
% colormap summer
% hold on
% plot3(xgdes, ygdes, zgdes+0.01, 'b')% A mountain brook!
% hold off




%...............................................%
%.................. FUNCTIONS ..................%
%...............................................%

function derivative = MyNPderiv(f, a, b, i)
    % number of subintervals
    N = length(x)-1;
    % preallocates vector to store derivative
    dy = zeros(size(x));
    % approximates derivative at lower bound using forward difference
    dy(1) = (y(2)-y(1))/(x(2)-x(1));
    % approximates derivative at upper bound using backward difference
    dy(N+1) = (y(N+1)-y(N))/(x(N+1)-x(N));
    % approximates derivatives at all other nodes using central differences
    for i = 2:N
        dy(i) = (y(i+1)-y(i-1))/(x(i+1)-x(i-1));
    end
    
end


function [dy,x] = derivative(x,y,x_star)
    % number of subintervals
    N = length(x)-1;
    % preallocates vector to store derivative
    dy = zeros(size(x));
    % approximates derivative at lower bound using forward difference
    dy(1) = (y(2)-y(1))/(x(2)-x(1));
    % approximates derivative at upper bound using backward difference
    dy(N+1) = (y(N+1)-y(N))/(x(N+1)-x(N));
    % approximates derivatives at all other nodes using central differences
    for i = 2:N
        dy(i) = (y(i+1)-y(i-1))/(x(i+1)-x(i-1));
    end
    
    % approximates derivative at specified points via linear interpolation
        
end







% This function (numerically) calculates an approximation
% of the partial derivative of the funtion f(x,y) at the
% point (a,b).

% ------------------INPUT-------------------%
% A function f to differentiate
% A value a for the x-variable.
% A vlaue b for the y-variable.
% A variable i to choose which partial derivative to
% calculate.
% ------------------OUTPUT------------------%
% The value of the partial derivative.



function [x,y] = MyNewton(f, g, a, b)
% This function finds a solution (x,y) to the system of
% equations f(x,y) = 0 and g(x,y) = 0 using Newton's method
% for functions of several variables. 

% ------------------INPUT-------------------%
% Two functions f and g.
% An initial guess (a,b) for the solution (x,y).
% ------------------OUTPUT------------------%
% x and y coordinates for the found solution.

format long 
f = f(x,y) = 0;
g = g(x,y) = 0;

z = f(a,b)+





end




function [x,y] = MyCritical(f, a, b)
% This function finds a critical point of the function f(x,y).

% ------------------INPUT-------------------%
% A function f to maximize.
% An initial geuss (a,b) for the critical point (x,y).
% ------------------OUTPUT------------------%
% x and y coordinates for the found critical point.





end




%...............................................%
%........... The topography function ...........%
%...............................................%

% Do not make any changes to the following function!

% "Topgraphy" describes the topography, i.e. the terain,
% of a mountain. You will be searching for peaks and
% valleys on this mountain with different numerical and
% graphical methods.

function z = topography(x,y)
% This function calculates a height z of a point on a mountain
% terain, provided the coordinates x and y.

% This function works for matricies as x and y, and is thus
% compatible with meshgrid.

% ------------------INPUT-------------------%
% Value(s) for x and y.
% ------------------OUTPUT------------------%
% Value(s) for z.

z= 0.001*(y-4).^2   +0.001*(x-3).^2 ...   
+0.5*exp(-0.25.*(y-11).^2-0.01*x.^2)...     
+0.5*exp(-0.01.*y.^2-0.25*(x-11).^2)...     
+0.5*exp(-0.1.*(y+1).^2-0.01*(x-5).^2)...    
+0.3*exp(-0.1.*(y-3.5).^2-0.5*(x+1).^2)... 
+0.5*exp(-0.1.*(y-8).^2-0.1*(x-0).^2)...  
+1.*exp(-0.1.*(y-9).^2-0.1*(x-8.5).^2)...  
+0.5*exp(-0.25.*(y-6).^2-0.1*(x-6).^2)...   
+0.25*exp(-0.5.*(y-3).^2-0.5*(x-8).^2)...  
+0.5*exp(-(y-5).^2-0.5*(x-5).^2)...    
+0.25*exp(-0.75.*(y-2).^2-(x-8).^2)...
+0.5*exp(-(y-6).^2-0.5*(x-3).^2)...
+0.5*exp(-(y-5).^2-0.5*(x-9).^2)...
+0.5*exp(-(y-9).^2-0.5*(x-5).^2);

end
