%%
clc
clear all

%% Question 1
% Define the function that we need to optimize in terms of x and y
syms x y
z = ((1.5* x - 2.5)^2) + ((y + 2)^2) - x*y;


% The initial values that will be used for univariate search
x_i = 0; 
y_i = 5;

% Defining a function called univariate search that we call for finding the
% minimum of the function z
function [x_min, y_min] = univariate_search(z, x_i, y_i)
iter = 0;
max_iter = 1000;
x_old = x_i;
y_old = y_i;
err_tolerance = .5 * (10^(-1));
error = 1000;
% Defining a counter that decides on which direction to optimize:
% the method optimizes across x, when the x_or_y_counter is an even number,
% and optimizes across y, when the x_or_y_counter is an odd number, which
% makes it oscillate between the optimization across x and y
x_or_y_counter = 0;  

syms x y
z = z;
% Defining the ranges of plotting for x and y
x_plot = linspace(-3,7);
y_plot = linspace(-3,7);
[X,Y] = meshgrid(x_plot,y_plot);
Z = ((1.5* X - 2.5).^2) + ((Y + 2).^2) - X.*Y;
hold on; 
% Starting plotting the contour maps of the function z
contour(X,Y,Z,20)
% Plotting the initial point in black
plot(x_old, y_old, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k');

while iter < max_iter && error > err_tolerance
if rem(x_or_y_counter, 2) == 0

    z_x = matlabFunction(subs(z, y, y_old), 'Vars', x);  % This creates a numeric function of x only
    
    % Using fminsearch to minimize z_x with respect to x, doing the
    % x-dimentional optimization
    x_new = fminsearch(z_x, x_i);
    y_new = y_old;
    error = abs((subs(z,{x,y}, {x_new, y_new}) - subs(z,{x,y}, {x_old, y_old}))/ subs(z,{x,y}, {x_new, y_new}))*100;
    x_old = x_new;
    y_old = y_new;
     % Incrementing the counter
    x_or_y_counter = x_or_y_counter+1;
     % Plotting the points along the x-direction in red
    plot(x_new, y_new, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');  
else
    z_y = matlabFunction(subs(z, x, x_old), 'Vars', y);  % This creates a numeric function of x only
    
    % Using fminsearch to minimize z_y with respect to y, doing the
    % y-dimentional optimization
    y_new = fminsearch(z_y, y_i);
    x_new = x_old;
    error = abs((subs(z,{x,y}, {x_new, y_new}) - subs(z,{x,y}, {x_old, y_old}))/ subs(z,{x,y}, {x_new, y_new}))*100;
    x_old = x_new;
    y_old = y_new;
    % Incrementing the counter
    x_or_y_counter = x_or_y_counter+1;
    % Plotting the points along the y-direction in blue
    plot(x_new, y_new, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b');

end

end
x_min = x_new;
y_min = y_new;
end
[x_min, y_min] = univariate_search(z, x_i, y_i)
f_minimum = double(subs(z,{x,y}, {x_min, y_min}))
% Plotting the minimum point in green
plot(x_min, y_min, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'g');
xlabel('X')
ylabel('Y')
hold off;
%% Questtion 2

% Defining Newton's method that will be used to optimize the step h, the
% step we walk along the gradient direction
function h_new = newton_optimization(func, h_old)
    % Defining the error tolerance value
    e_s = .5 * (10^(2 - 3));  
    % Assining a big value for the initial error just for the first
    % iteration to run
    err = 1000;  
    iter = 0; 
    % Maximum iterations to make sure the method will not run infinitely
    max_iter = 1000;  
    syms h;
    % Passing the function that we need its maximum to f
    f = func;  
    % Calculating the first derivative of the function
    f_1 = diff(f, h);  
    % Calculating the second derivative of the function
    f_2 = diff(f_1, h);  
    % Converting the first and the second derivatives into matlab function
    % to easily handle substitution in them
    f_1_f = matlabFunction(f_1, 'Vars', h);
    f_2_f = matlabFunction(f_2, 'Vars', h);

    % Initialize x_new with the initial guess
    h_new = h_old; 

    % Calculate the new x value using the Newton-Raphson method, using the
    % search for values that make the first derivative equal to zero (extreme value)
    while err > e_s && iter < max_iter
        h_new = h_old - (f_1_f(h_old) / f_2_f(h_old));  
        % Calculating the error
        err = abs((h_new - h_old) / h_new) * 100;
        % Updating the h_old value for the next iteration
        h_old = h_new;  
        iter = iter + 1; 
    end
end
% Defining the function that needs to be optimized
syms x y;
f = ((1.5 * x - 2.5)^2) + ((y + 2)^2) - x * y;

% Define the starting point for gradient descent
x_i = 0;
y_i = 5;


% Defining the ranges of plotting for x and y
x_plot = linspace(-3,7);
y_plot = linspace(-3,7);
[X,Y] = meshgrid(x_plot,y_plot);
Z = ((1.5* X - 2.5).^2) + ((Y + 2).^2) - X.*Y;
hold on; 
% Starting plotting the contour maps of the function z
contour(X,Y,Z,20)
% Plotting the initial point in black
plot(x_i, y_i, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k');

% Define the gradient descent method function
function [x_min, y_min] = steepest_descent(f, x_i, y_i)
iter = 0;
max_iter = 1000;
x_old = x_i;
y_old = y_i;
err_tolerance = .5 * (10^(2 - 3));
error = 1000;
while iter < max_iter && error > err_tolerance
    % Calculate the gradients with respect to x and y
    syms x y;
    f = f;
    grad_f_x = diff(f, x);
    grad_f_y = diff(f, y);
    % Evaluating the gradients at the point (x_old, y_old)
    grad_x_value = subs(grad_f_x, {x, y}, {x_old, y_old});
    grad_y_value = subs(grad_f_y, {x, y}, {x_old, y_old});
    
    % Updating the variables in terms of h multiplied by gradients of the
    % function with respect to x and y
    syms h;
    x_new = x_old - grad_x_value * h;
    y_new = y_old - grad_y_value * h;
    
    % Defining the function that will be optimized in terms of h (step size along the gradient direction) instead of
    % x and y
    g = ((1.5 * x_new - 2.5)^2) + ((y_new + 2)^2) - x_new * y_new;
    % Minimizing h using Newton's Method
    h_min = newton_optimization(g, 0);
    % Rcalculating the x, y values after h was minimized
    x_new = x_old - grad_x_value * h_min;
    y_new = y_old - grad_y_value * h_min;

    % Plotting the new point in red
    plot(x_new, y_new, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');  

    error = abs((subs(f,{x,y}, {x_new, y_new}) - subs(f,{x,y}, {x_old, y_old}))/ subs(f,{x,y}, {x_new, y_new}))*100;
    x_old = x_new;
    y_old = y_new;
    iter = iter + 1;
end
x_min = double(x_new);
y_min = double(y_new);
% Plotting the minimum point in green
plot(x_min, y_min, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'g');
xlabel('X')
ylabel('Y')
title('Contour Maps with the steps marked in red');
hold off;
end 

[x_min, y_min] =  steepest_descent(f, x_i, y_i)
double(subs(f, {x, y}, {x_min, y_min}))


% Create a mesh grid for x and y values for plotting
[x_grid, y_grid] = meshgrid(-5:0.1:5, -5:0.1:5);

% Evaluate the function on the mesh grid
f_grid = double(subs(f, {x, y}, {x_grid, y_grid}));
% Plot the surface
figure;
surf(x_grid, y_grid, f_grid, 'EdgeColor', 'none');

% Label the axes
colorbar;
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
title('Surface of the Function');

% Mark the minimum on the surface plot
hold on;
plot3(x_min, y_min, double(subs(f, {x, y}, {x_min, y_min})), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
hold off;


%% Question 3: Part b

% Define the range for x_1 and x_2 to create a plot
x_1 = linspace(0, 1000, 5000); 
x_2 = linspace(0, 1000, 5000); 
[X_1, X_2] = meshgrid(x_1, x_2); 

% Define the first constraint and giving z1 value 1 in case the constraint
% is fullfilled
% The size of the array z1 should be equal to the number of points plotted,
% X
Z1 = NaN(size(X_2)); 
Z1(5*X_1 + 20*X_2 <= 6000) = 1; 

% Define the second constraint and giving z2 value 1 in case the constraint
% is fullfilled
% The size of the array z2 should be equal to the number of points plotted,
% X
Z2 = NaN(size(X_2)); 
Z2(.14*X_1 + .05*X_2 <= 60) = 1; 

% Define the third constraint and giving z3 value 1 in case the constraint
% is fullfilled
% The size of the array z3 should be equal to the number of points plotted,
% X
Z3 = NaN(size(X_2)); 
Z3(2*X_1 + 25*X_2 <= 7000) = 1; 

% Create the intersection array that will be plotted as intersection area
% on the graph
Z_intersection = Z1 + Z2 + Z3; 
Z_intersection(Z_intersection < 3) = 1; 

% Creating the plot
figure;
hold on;

% Plotting the filled area for the first Constraint
surf(X_1, X_2, Z1, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Plotting the filled area for the second Constraint
surf(X_1, X_2, Z2, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Plotting the filled area for the third Constraint
surf(X_1, X_2, Z3, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Plotting the filled area for the intersection
surf(X_1, X_2, Z_intersection, 'FaceColor', 'y', 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Yellow for intersection

% Set the view and limits
%view(2); % View from the top
xlim([0 1000]); 
% Limiting the area we can see from x and y because 5000 points is actually
% so much to view
ylim([0 1000]); 
xlabel('Shovels');
ylabel('Blowers');
title('Regions satisfying the constraints and their intersection, with the optimum circled');
grid on;

% Calculating and plotting the boundary lines for the Constraints
boundary_x1 = linspace(0, 1200, 100); 
% Rearranging the constraint to obtain x_2 in terms of x_1
boundary_y1 = (6000 - 5*boundary_x1) / 20; 
% Keep y values within limits by disregarding y<0 or y>1000
boundary_y1(boundary_y1 < 0) = NaN; 
boundary_y1(boundary_y1 > 1000) = NaN; 
plot(boundary_x1, boundary_y1, 'r', 'LineWidth', 2); 

% Calculating and plot the boundary lines for the Constraints
boundary_x2 = linspace(0, 1200, 100); 
% Rearranging the constraint to obtain x_2 in terms of x_1
boundary_y2 = (60 - .14*boundary_x2) / .05; 
% Keep y values within limits by disregarding y<0 or y>1000
boundary_y2(boundary_y2 < 0) = NaN; 
boundary_y2(boundary_y2 > 1000) = NaN; 
plot(boundary_x2, boundary_y2, 'm', 'LineWidth', 2); 
% Calculating and plot the boundary lines for the Constraints
boundary_x3 = linspace(0, 1200, 100); 
% Rearranging the constraint to obtain x_2 in terms of x_1
boundary_y3 = (7000 - 2*boundary_x3) / 25; 

% Keep y values within limits by disregarding y<0 or y>1000
boundary_y3(boundary_y3 < 0) = NaN; 
boundary_y3(boundary_y3 > 1000) = NaN; 
plot(boundary_x3, boundary_y3, 'c', 'LineWidth', 2); 

% Calculate the intersection point of boundaries 2 and 3, since they are
% the outermost limits of the intersection region
syms x1 x2
eq1 = 0.14*x1 + 0.05*x2 == 60;
eq2 = 5*x1 + 20*x2 == 6000;
solution = vpasolve([eq1, eq2], [x1, x2]);
x_intersect = double(solution.x1);
y_intersect = double(solution.x2);

% Marking the intersection point on the graph with a black circle 
plot(x_intersect, y_intersect, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k'); 
fprintf('The number of shovels is %.2f (~%.2f)\n', x_intersect, round(x_intersect));
fprintf('The number of blowers is %.2f (~%.2f)\n', y_intersect, round(y_intersect));

hold off;

%% Question 3: Part c

% Solving the problem using matlab built-in linear programming
% Coefficients for constraints
A = [5, 20;    
     0.14, 0.05;
     2, 25;
     0, -1;    
     -1, 0]; 
% This matrix contains the other side of the inequalities
b = [6000;    
     60;
     7000;
     0;
     0];

% Using a negative sign for the coefficients of the objective function
% since linprog minimizes not maximizes
f = -[15; 45]; 
x = linprog(f, A, b)

fprintf('The number of shovels is %.2f (~%.2f)\n', x(1), round(x(1)));
fprintf('The number of blowers is %.2f (~%.2f)\n', x(2), round(x(2)));


%% Question 4: Easier Part
% Loading the data from the file to arrays representing the independent
% variables and the output
data = load('C:\Users\DELL\Documents\MATLAB\HW_8_dataset_comp_meth.mat');
temp = data.temp;
power_output = data.power_output;
n = length(power_output);
% Defining the threthold error so that s_r is correct to 6 significant
% figures
err_trethold = .5 * (10^(2 - 6));

% Defining the residuals and the s_r  function in terms of x=a_0 and y=a_1
syms x y;
residuals = power_output - (x + y * temp);
s_r = sum(residuals.^2); 

% Defining Newton's method that will be used to optimize the step h, the
% step we walk along the gradient direction
function h_new = newton_optimization_new(func, h_old)
    e_s = 0.5;  
    err = 100;  
    iter = 0; 
    max_iter = 1000;  
    syms h;
    
    f = func;  
    f_1 = diff(f, h);  
    f_2 = diff(f_1, h);  
    f_1_f = matlabFunction(f_1, 'Vars', h); 
    f_2_f = matlabFunction(f_2, 'Vars', h);  

    h_new = h_old; 

    while err > e_s && iter < max_iter
        h_new = h_old - (f_1_f(h_old) / f_2_f(h_old));  
        err = abs((h_new - h_old) / h_new) * 100; 
        h_old = h_new;  
        iter = iter + 1; 
    end
end

% Defining the function that needs to be optimized (minimized)
f = s_r;

% The initial guess for optimization we will use
x_i = 0;
y_i = 0;

% Defining the steepest descent method function, the same as in question 2
function [x_min, y_min] = steepest_descent_new(f, x_i, y_i)

    data = load('C:\Users\DELL\Documents\MATLAB\HW_8_dataset_comp_meth.mat');
    temp = data.temp;
    power_output = data.power_output;
    n = length(power_output);
    err_trethold = .5 * (10^(-4));
    % initialising the error for the first iteration to run
    err = 1000; 

    iter_n = 0;
    max_iter_n = 1000;
    x_old = x_i;
    y_old = y_i;
    

    while iter_n < max_iter_n && err > err_trethold
        syms x y;
        grad_f_x = diff(f, x);
        grad_f_y = diff(f, y);
        
        grad_x_value = subs(grad_f_x, {x, y}, {x_old, y_old});
        grad_y_value = subs(grad_f_y, {x, y}, {x_old, y_old});
        
        syms h;
        x_new = x_old - grad_x_value * h;
        y_new = y_old - grad_y_value * h;
        
        residuals = power_output - (x_new + y_new * temp);
        % Defining a function g that will be optimized but in terms of the
        % variable h instead of x and y
        g = sum(residuals.^2);

        h_min = newton_optimization_new(g, 0);
        x_new = x_old - grad_x_value * h_min;
        y_new = y_old - grad_y_value * h_min;

        residuals_old = power_output - (x_old + y_old * temp);
        sr_old = sum(residuals_old.^2);
        residuals_new = power_output - (x_new + y_new * temp);
        sr_new = sum(residuals_new.^2);

        x_old = x_new;
        y_old = y_new;

        err = double(abs(((sr_new - sr_old)/sr_new))*100);
        iter_n = iter_n + 1;
    end
    x_min = double(x_new);
    y_min = double(y_new);
end 

[a_0, a_1] = steepest_descent_new(f, x_i, y_i)

% Plotting the scattered data points and the line we found to see 
hold on;
plot(data.temp, data.power_output, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');

% Defining of the straight line that best fit the data after finding a_0 and a_1
f_optm = a_0 + a_1 * data.temp;
plot(data.temp, f_optm, 'b')
xlabel('Temperature');
ylabel('Power Output');