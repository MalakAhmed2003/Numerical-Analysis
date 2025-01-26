%%
clc
clear all
%% Question 1: Part a
% Defining the function to be graphed
f = @(x) 4 - 9*x - 2*(x.^2) - 8*(x.^3) + 2*(x.^4);  
% Asigning the range of x that the we will graph the function on
x = linspace(3.2, 3.3, 10000);  

% Plotting the function 
plot(x, f(x), 'b'); 
xlabel('x');
ylabel('f(x)');
title('Plot of f(x)');
grid on;

%% Question 1: Part b

% Define the symbolic variable x that will be the independent variable of
% the function to be minimized
syms x; 
% Defining the function
f = 4 - 9*x - 2*(x^2) - 8*(x^3) + 2*(x^4);  
% Defining the nitial guess that we will use for the search of the maximum
x_0 = 1; 

% Defining the Newton optimization function
function x_new = newton_optimization(func, x_old)
    % Defining the error tolerance value
    e_s = 0.5;  
    % Assining a big value for the initial error just for the first
    % iteration to run
    err = 100;  
    iter = 0; 
    % Maximum iterations to make sure the method will not run infinitely
    max_iter = 1000;  
    syms x;
    % Passing the function that we need its maximum to f
    f = func;  
    % Calculating the first derivative of the function
    f_1 = diff(f, x);  
    % Calculating the second derivative of the function
    f_2 = diff(f_1, x);  
    % Converting the first and the second derivatives into matlab function
    % to easily handle substitution in them
    f_1_f = matlabFunction(f_1, 'Vars', x) 
    f_2_f = matlabFunction(f_2, 'Vars', x)  

    % Initialize x_new with the initial guess
    x_new = x_old; 

    % Calculate the new x value using the Newton-Raphson method, using the
    % search for values that make the first derivative equal to zero (extreme value)
    while err > e_s && iter < max_iter
        x_new = x_old - (f_1_f(x_old) / f_2_f(x_old))  % Newton's method definition
        % Calculating the error
        err = abs((x_new - x_old) / x_new) * 100 
        % Updating the x_old value for the next iteration
        x_old = x_new;  
        iter = iter + 1; 
    end
end
% Calling the Newton optimization function with the initial guess that we
% already defined
x_minimum = newton_optimization(f, x_0)
y = double(subs(f,x_minimum))
%% Question 2: Part a

% Defining the function to be graphed
f = @(x) (1./(.03+(x-.8).^2)) + (1./(.05+(x-2).^2))-7; 
% Asigning the range of x that the we will graph the function on
x = linspace(-10, 10, 10000);  

% Plotting the function
plot(x, f(x), 'b'); 
xlabel('x');
ylabel('f(x)');
title('Plot of f(x)');
grid on;
%% Question 2: Part b
% Defining the function we need its maximum
f = @(x) (1./(.03+(x-.8).^2)) + (1./(.05+(x-2).^2))-7;  
% Creating a matlab function that applies the golden-search algorithm
function [x_opt] = golden_seach(f, x_u, x_l)

iter_max = 1000;
% Defining the error tolerance value
e_s = .00001;
% Defining the golden ratio
R = .618;
 % Assining a big value for the initial error just for the first
 % iteration to run
err =1000;
iter =0;
d = (x_u - x_l) * R;
while err > e_s && iter < iter_max
    x_1 = x_l +d;
    x_2 = x_u -d;
    if x_1 > x_2
        if f(x_1) > f(x_2)
            x_opt = x_1;
            x_l = x_2;
            err = (1-R)*abs((x_u - x_l)*100/x_opt);

        end 
        if f(x_2) > f(x_1)
            x_opt = x_2;
            x_u = x_1;
            err = (1-R)*abs((x_u - x_l)*100/x_opt);
            
        end 
    end
    if x_2 > x_1
        if f(x_1) > f(x_2)
            x_opt = x_1;
            x_u = x_2;
            err = (1-R)*abs((x_u - x_l)*100/x_opt);

        end 
        if f(x_2) > f(x_1)
            x_opt = x_2;
            x_l = x_1;
            err = (1-R)*abs((x_u - x_l)*100/x_opt);
            
        end 
    end
    d = (x_u - x_l) * R;

    iter = iter +1;
end
end
 x_max_glob = golden_seach(f, 1.5, 0)
 y_max_glob = f(x_max_glob)
 %% Question 2: Part c

% Defining the function we need its maximum
f = @(x) (1./(.03+(x-.8).^2)) + (1./(.05+(x-2).^2))-7;  
% Creating a matlab function that applies the golden-search algorithm
function [x_opt] = golden_seach_local(f, x_u, x_l)

iter_max = 1000;
% Defining the error tolerance value
e_s = .00001;
% Defining the golden ratio
R = .618;
% Assining a big value for the initial error just for the first
% iteration to run
err =1000;
iter =0;
d = (x_u - x_l) * R;
while err > e_s && iter < iter_max
    x_1 = x_l +d;
    x_2 = x_u -d;
    if x_1 > x_2
        if f(x_1) > f(x_2)
            x_opt = x_1;
            x_l = x_2;
            err = (1-R)*abs((x_u - x_l)*100/x_opt);

        end 
        if f(x_2) > f(x_1)
            x_opt = x_2;
            x_u = x_1;
            err = (1-R)*abs((x_u - x_l)*100/x_opt);
            
        end 
    end
    if x_2 > x_1
        if f(x_1) > f(x_2)
            x_opt = x_1;
            x_u = x_2;
            err = (1-R)*abs((x_u - x_l)*100/x_opt);

        end 
        if f(x_2) > f(x_1)
            x_opt = x_2;
            x_l = x_1;
            err = (1-R)*abs((x_u - x_l)*100/x_opt);
            
        end 
    end
    d = (x_u - x_l) * R;

    iter = iter +1;
end
end
% Modifying the starting boundaries so that I can get the local maximum
% inst
 x_max_loca = golden_seach_local(f, 3, 1)
 y_max_loca = f(x_max_loca)
 %% Question 2: Part d

 % Defining the function to be graphed
f = @(x) (1./(.03+(x-.8).^2)) + (1./(.05+(x-2).^2))-7;  
% Asigning the range of x that the we will graph the function on
x = linspace(-10, 10, 10000);  

% Plotting the function itself
plot(x, f(x), 'b');  
hold on;
% Plotting the local and global maxima
plot(x_max_loca, y_max_loca, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');  
plot(x_max_glob, y_max_glob, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');  
text(x_max_loca, y_max_loca, sprintf(' Local Maximum(%.1f, %.1f)', ...
    x_max_loca, y_max_loca),'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
text(x_max_glob, y_max_glob, sprintf(' Global Maximum (%.1f, %.1f)', ...
    x_max_glob, y_max_glob),'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
xlabel('x');
ylabel('f(x)');
title('Plot of f(x)');
grid on;
%% Question 3:

%Defining the function symbolically
syms x y
f = 2*x + 3*y - 13*(x^2) - 7*(y^2) + x*y; % Define the function

% Convert the symbolic function to a matlab function to easily substitute
% and handle it
f_numeric = matlabFunction(f);
% Initializng f_max with a very small value for the first iteration to run
f_max = -1 * (10^9);
for i = 1:1000000
    % Guessing a random x in [-1.5, 1.5]
    r_x = -1.5 + 3 * rand(1, 1); 
    % Guessing a random y in [-3, 3]
    r_y = -3 + 6 * rand(1, 1);   
    
    % Evaluate the function with the random values we have just guessed
    if f_numeric(r_x, r_y) > f_max
        % Updating f_max only if a bigger value is found
        f_max = f_numeric(r_x, r_y); 
        x_max = r_x;
        y_max = r_y;
    end
end

sprintf('Maximum value found: %.1f', f_max)

% Create meshgrid for the surface variables x and y
[x_grid, y_grid] = meshgrid(linspace(-1.5, 1.5, 1000), linspace(-3, 3, 1000));
z = 2*x_grid + 3*y_grid - 13*(x_grid.^2) - 7*(y_grid.^2) + x_grid.*y_grid; 

% Plotting the surface
surf(x_grid, y_grid, z, 'EdgeColor', 'none');
hold on;

% Plot the maximum point
plot3(x_max, y_max, f_numeric(x_max, y_max), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(x_max, y_max, f_numeric(x_max, y_max), ...
    sprintf(' Global maximum (%.2f, %.2f, %.2f)', x_max, y_max, f_numeric(x_max, y_max)), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'k');
colorbar;
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Surface Plot of f(x,y)');
hold off;
%% Showing the guesses that were made

% Defining the function symbolically
syms x y
f = 2*x + 3*y - 13*(x^2) - 7*(y^2) + x*y; % Define the function

% Convert the symbolic function to a MATLAB function for easier handling
f_numeric = matlabFunction(f);

% Initializing f_max with a very small value for the first iteration
f_max = -1 * (10^9);


% Create meshgrid for the surface variables x and y
[x_grid, y_grid] = meshgrid(linspace(-1.5, 1.5, 1000), linspace(-3, 3, 1000));
z = f_numeric(x_grid, y_grid);  % Calculate z values directly using f_numeric

% Plotting the surface
surf(x_grid, y_grid, z, 'EdgeColor', 'none');
hold on;

for i = 1:1000
    % Guessing a random x in [-1.5, 1.5]
    r_x = -1.5 + 3 * rand(1, 1); 
    % Guessing a random y in [-3, 3]
    r_y = -3 + 6 * rand(1, 1);   
    
    % Evaluate the function with the random values we have just guessed
    current_value = f_numeric(r_x, r_y);
    % Plot random points in black
    plot3(r_x, r_y, current_value, 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

    if current_value > f_max
        % Updating f_max only if a bigger value is found
        f_max = current_value; 
        x_max = r_x;
        y_max = r_y;
    end
end

sprintf('Maximum value found: %.1f', f_max);

% Plot the maximum point in red
plot3(x_max, y_max, f_numeric(x_max, y_max), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
text(x_max, y_max, f_numeric(x_max, y_max) + 0.5, ...  
    sprintf(' Global maximum (%.2f, %.2f, %.2f)', x_max, y_max, f_numeric(x_max, y_max)), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'k');

colorbar;
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Surface Plot of f(x,y)');
hold off;
