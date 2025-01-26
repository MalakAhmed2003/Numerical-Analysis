%%
clc
clear all
%% Question 1: Part 1

% Defining the height function that will be optimized
syms x1 x2
f_x1_x2 = (1000 /((abs(x1 - 1.6) + abs(x2 - 1.6) + 1))) ...
+ ( 2000 / ((abs(x1 - 2.1) ...
+ abs(x2 - 2.1) + 1))) +(2000 / ((abs(x1 - 1.1) + abs(x2 + 3.1) + 1))) ...
+ (4000 / (abs(x1 - 0.1) + abs(x2 - 0.1) + 1));

% Converting the symbolic function to a matlab function to easily handle it
f_numeric = matlabFunction(f_x1_x2, 'Vars', {x1, x2});

% Defining the range for x1 and x2 (Chosen based on the sytematic way the
% critical points tend to appear near values that make f' undefined)
x1_vals = linspace(-3, 3, 1000);  
x2_vals = linspace(-4, 4, 1000);  

% Creating a meshgrid for x1 and x2
[X1, X2] = meshgrid(x1_vals, x2_vals);

% Compute the values of the function at each point on the grid
Z = f_numeric(X1, X2);

% Creating a surface of the function
figure;
surf(X1, X2, Z);
shading interp;  
xlabel('x1');
ylabel('x2');
zlabel('f(x1, x2)');
title('Surface plot of f(x1, x2)');
colorbar;
grid on;

%% Question 1: Part 2

% Part a
% Defining the -ve height function that will be optimized
f_x1_x2 = @(x1, x2) (-1000 / (abs(x1 - 1.6) + abs(x2 - 1.6) + 1)) + ...
                    (-2000 / (abs(x1 - 2.1) + abs(x2 - 2.1) + 1)) + ...
                    (-2000 / (abs(x1 - 1.1) + abs(x2 + 3.1) + 1)) + ...
                    (-4000 / (abs(x1 - 0.1) + abs(x2 - 0.1) + 1));

% Initial guess for optimization (to find global max)
Initial_guess_1 = [0, 0];

% Usig fminsearch to find the minimum (which corresponds to the global maximum in the original function)
[optimal_guess_1, f_min_1] = fminsearch(@(x) f_x1_x2(x(1), x(2)), Initial_guess_1);
f_min_1 = f_min_1/-1;

disp('The global maximum coordinates are:');
disp(optimal_guess_1);
disp('Function Value at the global maximum is:');
disp(f_min_1);

% Initial guess for optimization (local max)
Initial_guess_2 = [1.4, 1.4];

% Using fminsearch to find the minimum (which corresponds to a local maximum in the original function)
[optimal_guess_2, f_min_2] = fminsearch(@(x) f_x1_x2(x(1), x(2)), Initial_guess_2);
f_min_2 = f_min_2/-1;


disp('The first local maximum coordinates are:');
disp(optimal_guess_2);
disp('Function Value at the first local maximum is:');
disp(f_min_2);

% Initial guess for optimization (local max)
Initial_guess_3 = [1, -3];

% Using fminsearch to find the minimum (which corresponds to a local maximum in the original function)
[optimal_guess_3, f_min_3] = fminsearch(@(x) f_x1_x2(x(1), x(2)), Initial_guess_3);
f_min_3 = f_min_3/-1;

disp('The second local maximum coordinates are:');
disp(optimal_guess_3);
disp('Function Value at the second local maximum is:');
disp(f_min_3);

% Initial guess for optimization (local max)
Initial_guess_4 = [1.9, 1.9];

% Using fminsearch to find the minimum (which corresponds to a local maximum in the original function)
[optimal_guess_4, f_min_4] = fminsearch(@(x) f_x1_x2(x(1), x(2)), Initial_guess_4);
f_min_4 = f_min_4/-1;

disp('The third local maximum coordinates are:');
disp(optimal_guess_4);
disp('Function Value at the third local maximum is:');
disp(f_min_4);

% Part b:

syms x1 x2
f_symbolic = (1000 /((abs(x1 - 1.6) + abs(x2 - 1.6) + 1))) + ...
              (2000 / ((abs(x1 - 2.1) + abs(x2 - 2.1) + 1))) + ...
              (2000 / ((abs(x1 - 1.1) + abs(x2 + 3.1) + 1))) + ...
              (4000 / (abs(x1 - 0.1) + abs(x2 - 0.1) + 1));


% Plotting Procedures just like part 1
f_numeric = matlabFunction(f_symbolic, 'Vars', {x1, x2});


x1_vals = linspace(-3, 3, 500);  % Range for x1
x2_vals = linspace(-4, 4, 500);  % Range for x2


[X1, X2] = meshgrid(x1_vals, x2_vals);


Z = f_numeric(X1, X2);


figure;
surf(X1, X2, Z);
shading interp;  
xlabel('x1');
ylabel('x2');
zlabel('f(x1, x2)');
title('Surface plot of f(x1, x2)');
colorbar;
grid on;

hold on;
% Marking out the global maximum
plot3(optimal_guess_1(1), optimal_guess_1(2), f_min_1, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(optimal_guess_1(1), optimal_guess_1(2), f_min_1, ...
    sprintf('Global maximum (%.2f, %.2f, %.2f)', optimal_guess_1(1), optimal_guess_1(2), f_min_1), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'k');

% Marking out the other local maximum
plot3(optimal_guess_2(1), optimal_guess_2(2), f_min_2, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(optimal_guess_2(1), optimal_guess_2(2), f_min_2, ...
    sprintf('Local maximum (%.2f, %.2f, %.2f)', optimal_guess_2(1), optimal_guess_2(2), f_min_2), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'k');

% Marking the second local maximum
plot3(optimal_guess_3(1), optimal_guess_3(2), f_min_3, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(optimal_guess_3(1), optimal_guess_3(2), f_min_3, ...
    sprintf('Local maximum (%.2f, %.2f, %.2f)', optimal_guess_3(1), optimal_guess_3(2), f_min_3), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'k');

% Marking the third local maximum
plot3(optimal_guess_4(1), optimal_guess_4(2), f_min_4, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(optimal_guess_4(1), optimal_guess_4(2), f_min_4, ...
    sprintf('Local maximum (%.2f, %.2f, %.2f)', optimal_guess_4(1), optimal_guess_4(2), f_min_4), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'k');



%% Question 1: Part 3 - Plotting Surface with Constraints 

% Part a
% Defining the -ve height function that will be optimized
f_x1_x2 = @(x1, x2) (-1000 / (abs(x1 - 1.6) + abs(x2 - 1.6) + 1)) + ...
                    (-2000 / (abs(x1 - 2.1) + abs(x2 - 2.1) + 1)) + ...
                    (-2000 / (abs(x1 - 1.1) + abs(x2 + 3.1) + 1)) + ...
                    (-4000 / (abs(x1 - 0.1) + abs(x2 - 0.1) + 1));

% Defining the matrices of the constrains system: A * x <= b
A = [-1.2 2;
      0.6 1];
b = [-2;
      1.5];

% Initial guess for optimization
x0 = [1, -3];  

% Calling fmincon to find the minimum (which corresponds to the maximum in the original function)
x_opt = fmincon(@(x) f_x1_x2(x(1), x(2)), x0, A, b);

% Part b
y_opt = f_x1_x2(x_opt(1),x_opt(2) )/-1;

disp('Optimal constrained solution (x1, x2) is:');
disp(x_opt);
disp('The value of f at the optimal solution is:');
disp(y_opt);


% Part c
% Plotting procedures as in part 1
syms x1 x2
f_symbolic = (1000 / (abs(x1 - 1.6) + abs(x2 - 1.6) + 1)) + ...
              (2000 / (abs(x1 - 2.1) + abs(x2 - 2.1) + 1)) + ...
              (2000 / (abs(x1 - 1.1) + abs(x2 + 3.1) + 1)) + ...
              (4000 / (abs(x1 - 0.1) + abs(x2 - 0.1) + 1));

% Converting the symbolic function to a matlab numerical function
f_numeric = matlabFunction(f_symbolic, 'Vars', {x1, x2});

% Defining the range for x1 and x2
x1_vals = linspace(-3, 3, 500);  
x2_vals = linspace(-4, 4, 500);  

% Creating a meshgrid for x1 and x2
[X1, X2] = meshgrid(x1_vals, x2_vals);

% Computing the values of the function at each point on the grid
Z = f_numeric(X1, X2);

% Creating a surface plot of the function
figure;
surf(X1, X2, Z);
shading interp;  
xlabel('x1');
ylabel('x2');
zlabel('f(x1, x2)');
title('Surface plot of f(x1, x2)');
colorbar;
grid on;
hold on;

% Defining the range for x and z for the constraints plotting
x_const_1 = linspace(-3, 3, 500);   
z_const_1 = 1000:100:6000;  

x_const_2 = linspace(-3, 3, 500);   
z_const_2 = 1000:100:6000;   

% Creating meshgrids for the constraint surfaces
[X_const_1, Z_const_1] = meshgrid(x_const_1, z_const_1);
[X_const_2, Z_const_2] = meshgrid(x_const_2, z_const_2);

% Computing corresponding y values for each constraint surface
Y_const_1 = (1.2 * X_const_1 - 2) / 2;  
Y_const_2 = (-0.6 * X_const_1 + 1.5);   

% Plottig the constraint surfaces
surf(X_const_1, Y_const_1, Z_const_1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');  
surf(X_const_2, Y_const_2, Z_const_2, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 


% Plotting the optimal solution found by fmincon within the constrained
% area
plot3(x_opt(1), x_opt(2), f_numeric(x_opt(1), x_opt(2)), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(x_opt(1), x_opt(2), f_numeric(x_opt(1), x_opt(2)), ...
    sprintf('  Optimal (%.2f, %.2f)', x_opt(1), x_opt(2)), ...
    'FontSize', 10, 'Color', 'k');


legend({'Objective Function Surface', 'Constraint 1', 'Constraint 2', 'Optimal Point'}, 'Location', 'best');
grid on;


%% Question 1: Part 4
% Part a

% Defining the function f_x1_x2
f_x1_x2 = @(x1, x2) (1000 / (abs(x1 - 1.6) + abs(x2 - 1.6) + 1)) + ...
                    (2000 / (abs(x1 - 2.1) + abs(x2 - 2.1) + 1)) + ...
                    (2000 / (abs(x1 - 1.1) + abs(x2 + 3.1) + 1)) + ...
                    (4000 / (abs(x1 - 0.1) + abs(x2 - 0.1) + 1));

% Defining the function g
g = @(x1, x2) (300 * (sqrt((x1 - 4.5)^2 + (x2 - 1.2)^2))^1.7);

% Initialize the number of sampled values (rows of the matrix)
n = 10000;
% Initialising the matrix that contains x1, x2, f (x1, x2) and g(x1, x2)
rand_values = zeros(n, 4);

% Set the initial values that we start the hike from
rand_values(1, 1) = 4.5;
rand_values(1, 2) = 1.2;
rand_values(1, 3) = f_x1_x2(4.5, 1.2);  
rand_values(1, 4) = g(4.5, 1.2);       

% Looping to fill in the matrix
for i = 2:n  
    % x values between -9 and 9 (the range where the function of the height has its maximums, so that we do not have to plot sooooo many points away from the peaks )
    rand_values(i, 1) = -9 + 18 * rand(1,1);
    % y values between -9 and 9
    rand_values(i, 2) = -9 + 18 * rand(1,1);
    % Evaluate the functions at the random values
    rand_values(i, 3) = f_x1_x2(rand_values(i, 1), rand_values(i, 2));  % f_x1_x2 value
    rand_values(i, 4) = g(rand_values(i, 1), rand_values(i, 2));        % g value
end

% Display the matrix with all the values
disp('The Matrix is')
disp(rand_values);

% Part b
x = double(rand_values(:, 4));  
y = double(rand_values(:, 3)); 

% Plotting the Pareto Plot
figure;
plot(x, y, 'LineStyle', 'none', 'Marker', '.', 'Color', [0 0 1], 'MarkerSize', 10);
xlabel('g(x1, x2)');
ylabel('f(x1, x2)');
title('The Pareto Plot');

% Part c

% Defining the prtp function to find the Pareto front
function [front] = prtp(x, y, n)
    % Indicating the initial front matrix
    front = []; 
    % Index for the Pareto front
    k = 1;  

    for i = 1:n
        include = true;  
        for j = 1:n
            % If another point dominates the current point, exclude it
            if (x(j) <= x(i)) && (y(j) >= y(i)) && (i ~= j)
                include = false;
                break;  
            end
        end
        
        % If the point is not dominated, we add it to the Pareto front
        % matrix
        if include
            % Storing the value of g
            front(k, 1) = x(i);
            % Storing the value of f
            front(k, 2) = y(i);  
            % Increasing k the index for the next Pareto point
            k = k + 1;  
        end
    end
end

A = prtp(x, y, n);

% Plotting the Pareto front of A as red dots
figure;
plot(x, y, 'LineStyle', 'none', 'Marker', '.', 'Color', [0 0 1], 'MarkerSize', 10);
hold on;
xlabel('g(x1, x2)');
ylabel('f(x1, x2)');
title('The Pareto Plot');
plot(A(:, 1), A(:, 2), 'LineStyle', 'none', 'Marker', '.', 'Color', [1 0 0], 'MarkerSize', 15);
hold off; 

%% Question 1: Part 5

% Part a

% Defining the weights for the objective functions
wt = [2/3, 1/3];  


syms x1 x2
% Defining the function t as a weighted sum of f and g
t = wt(1)*((1000 / (abs(x1 - 1.6) + abs(x2 - 1.6) + 1)) + ...
           (2000 / (abs(x1 - 2.1) + abs(x2 - 2.1) + 1)) + ...
           (2000 / (abs(x1 - 1.1) + abs(x2 + 3.1) + 1)) + ...
           (4000 / (abs(x1 - 0.1) + abs(x2 - 0.1) + 1))) - ...
    wt(2)*(300 * (sqrt((x1 - 4.5)^2 + (x2 - 1.2)^2))^1.7);

% Convert the symbolic function to a matlab numerical function to easily
% handle it
f_numeric = matlabFunction(t, 'Vars', {x1, x2}); 

% Defining the ranges for x1 and x2
x1_vals = linspace(-5, 5, 500);  % Range for x1
x2_vals = linspace(-5, 5, 500);  % Range for x2

% Creating a meshgrid for x1 and x2
[X1, X2] = meshgrid(x1_vals, x2_vals);

Z = f_numeric(X1, X2);

% Creating a plot of the function t
figure;
surf(X1, X2, Z);  
hold on;
shading interp;  
xlabel('x1');
ylabel('x2');
zlabel('t(x1, x2)');
title('Surface plot of t(x1, x2) with the global max');
colorbar;
grid on;

% Initial guess for optimization
Initial_guess_1 = [.5, .5];

t_optim = @(x1, x2) (-1*(wt(1)*((1000 / (abs(x1 - 1.6) + abs(x2 - 1.6) + 1)) + ...
           (2000 / (abs(x1 - 2.1) + abs(x2 - 2.1) + 1)) + ...
           (2000 / (abs(x1 - 1.1) + abs(x2 + 3.1) + 1)) + ...
           (4000 / (abs(x1 - 0.1) + abs(x2 - 0.1) + 1))) - ...
    wt(2)*(300 * (sqrt((x1 - 4.5)^2 + (x2 - 1.2)^2))^1.7)));

% Using fminsearch to find the minimum (which corresponds to the global maximum in the original function)
[optimal_guess_1, f_min_1] = fminsearch(@(x) t_optim(x(1), x(2)), Initial_guess_1);
f_min_1 = f_min_1/-1;

disp('The global maximum coordinates of t(x1, x1) are:');
disp(optimal_guess_1);
disp('Function Value at that global maximum is:');
disp(f_min_1);
% Plotting the optimal point found
plot3(optimal_guess_1(1), optimal_guess_1(2), f_min_1, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(optimal_guess_1(1), optimal_guess_1(2), f_min_1, ...
    sprintf('Global maximum (%.2f, %.2f, %.2f)', optimal_guess_1(1), optimal_guess_1(2), f_min_1), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'k');
hold off;

% Part b
% Defining the matrices that will be used to solve the constrained system
A = [-1.2 2;
      0.6 1];
b = [-2;
      1.5];

% Initial guess for constrained optimization
x0 = [1, -3];  


% Calling fmincon to find the minimum (which corresponds to the maximum in the original function)
x_opt = fmincon(@(x) t_optim(x(1), x(2)), x0, A, b);
y_opt = t_optim(x_opt(1),x_opt(2) )/-1;

disp('Constrained optimal solution (x1, x2) is:');
disp(x_opt);
disp('The value of f at the constrained optimal solution is:');
disp(y_opt);



% Defining the range for x and z for the constraints
x_const_1 = linspace(-3, 3, 500);   
z_const_1 = 1000:100:6000;  

x_const_2 = linspace(-3, 3, 500);   
z_const_2 = 1000:100:6000;   

% Creating meshgrids for the constraint surfaces
[X_const_1, Z_const_1] = meshgrid(x_const_1, z_const_1);
[X_const_2, Z_const_2] = meshgrid(x_const_2, z_const_2);

% Computing corresponding y values for each constraint surface
Y_const_1 = (1.2 * X_const_1 - 2) / 2;  
Y_const_2 = (-0.6 * X_const_1 + 1.5);   

% Plotting the constraint surfaces 
figure;
surf(X1, X2, Z); 
hold on;
shading interp;  
xlabel('x1');
ylabel('x2');
zlabel('t(x1, x2)');
title('Surface plot of t(x1, x2) with the constrained optimum');
colorbar;
grid on;
surf(X_const_1, Y_const_1, Z_const_1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');  % Transparent surface
surf(X_const_2, Y_const_2, Z_const_2, 'FaceAlpha', 0.5, 'EdgeColor', 'none');  % Transparent surface


% Plot the constrained optimal solution found by fmincon on the constrained surface
plot3(x_opt(1), x_opt(2), f_numeric(x_opt(1), x_opt(2)), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(x_opt(1), x_opt(2), f_numeric(x_opt(1), x_opt(2)), ...
    sprintf('  Constrained Optimal (%.2f, %.2f)', x_opt(1), x_opt(2)), ...
    'FontSize', 10, 'Color', 'k');

legend({'Objective Function Surface', 'Constraint 1', 'Constraint 2', 'Optimal Point'}, 'Location', 'best');
grid on;


%% Question 2: Using trapizoidal method over the different intervals
X = [0 2 3 4 6 8 10];
p = [4.00 3.95 3.89 3.80 3.60 3.41 3.30];
A = (10^-4)*[105 108 109 113 125 137 150];
Y = p.*A;

n = length(p);
I =[];

for i= 1:1:n-1
    I(i) = (X(i+1)-X(i))*(Y(i)+Y(i+1))/2;
end
mass = 0;
for i= 1:1:n-1
    mass = mass +I(i);
end

mass = mass * (10^-3);
disp(mass)
%% Question 2: Using both Trapizoidal and Simpson's Method

X = [0 2 3 4 6 8 10];
p = [4.00 3.95 3.89 3.80 3.60 3.41 3.30];
A = (10^-4)*[105 108 109 113 125 137 150];
Y = p.*A;

I_new = [];
% Using Trpizoidal method for the first interval
for i= 1:1:1
    I(i) = (X(i+1)-X(i))*(Y(i)+Y(i+1))/2;
end
% Using Simpson's rule for three points, eveenly spaced from 2 to 4
I(2) = 2*(Y(2)+4*Y(3)+Y(4))/6;

% Using Simpson's rule for four points, evenly spaced from 4 to 10
I(3) = 6*(Y(4)+3*Y(5)+3*Y(6)+Y(7))/8;

mass = 0;
for i= 1:1:3
    mass = mass +I(i);
end

mass = mass * (10^-3);
disp('Mass is equal to (Kg):')
disp(mass)

%% Question 3: Part a

% Range of x
x = 0:0.1:0.4;   
% Number of steps
n = 5;
% Step size
h = 0.1;          
y_euler = zeros(1, n);  
z_euler = zeros(1, n);  

% Initial conditions
y_euler(1) = 2;         
z_euler(1) = 4;         

% Looping to update y and z using Euler's method
for i = 2:n
    % Calculate f_1 and f_2 at the current step
    f_1 = -y_euler(i-1) + 5*exp(x(i-1));        
    f_2 = (-y_euler(i-1) * z_euler(i-1)^2) / 4;       

    % Update y and z using Euler's method
    y_euler(i) = y_euler(i-1) + h * f_1;
    z_euler(i) = z_euler(i-1) + h * f_2;
end

%% Question 3: Part b

% Range of x
x = 0:0.1:0.4;   
n = 5; 
% Step size
h = 0.1;          
y = zeros(1, n);  
z = zeros(1, n);  

% Initial conditions
y(1) = 2;         
z(1) = 4;         

% Looping to update y and z using Midpoint method
for i = 1:n-1
    f_1_i = -y(i) + 5*exp(x(i));
    y_half = y(i) + f_1_i * h/2;

    % Calculating f_1 at the current step using the calculated y_half
    f_1 = -y_half + 5*exp(x(i)+(h/2));


    f_2_i = (-y(i) * z(i)^2) / 4; 
    z_half = z(i) + f_2_i * h/2;

     % Calculating f_1 at the current step using the calculated z_half
    f_2 = (-y_half * (z_half)^2) / 4;       

    % Updating y and z 
    y(i+1) = y(i) + h * f_1;
    z(i+1) = z(i) + h * f_2;
end

% Plot y vs x 
figure;  
yyaxis left
plot(x, y_euler, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'y using euler method');
hold on;
plot(x, y, '-x', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'y using midpoint method');
xlabel('x');
ylabel('y values');
grid on;

% Plot z vs x  
yyaxis right
plot(x, z_euler, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'z using euler method');
hold on;
plot(x, z, '-x', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'z using midpoint method');
xlabel('x');
ylabel('z values');
title('z and y vs x');
legend('show');
grid on;



%% Question 4
m = 1;         % mass (kg)
k = 9.8;       % spring constant (N/m^2)
g = 9.8;       % gravitational acceleration (m/s^2)
c = 0.5;       % damping coefficient (N s/m)

syms v(t) y(t) t

% Define the equations of motion
ode1 = diff(y) == v;  % y' = v (position is the integral of velocity)
ode2 = diff(v) == -g - (c/m)*v - (k/m)*y^2;  % v' = acceleration due to gravity, damping, and spring

odes = [ode1; ode2];
[VF,Sbs] = odeToVectorField(odes);
odsefcn = matlabFunction(VF, 'Vars',{t,y});
% Solve the system of equations
S = ode45(odesfcn);

% Display the solution
disp(S);



%% Question 4

% Initial conditions: [initial position, initial velocity]
x0 = [0.5, 0];  % Initial position is 0.5 m, initial velocity is 0 m/s

% Time vector: From t=0 to t=10 seconds, with 40*60 steps (i.e., 2400 time points)
t = linspace(0, 10, 40*60);  

% Solve the ODE using ode23s
[t, x] = ode23s(@OneDOFMSD_ex, t, x0);

% Plot the position versus time (x(:,1) is the position)
figure;
plot(t, x(:,1));  % Plot the first column of x, which is the position y(t)
xlabel('Time (s)');  % Label the x-axis
ylabel('Position (m)');  % Label the y-axis
title('Position vs. Time for the 1D Mass-Spring-Damper System');  % Plot title
grid on;  % Add grid to the plot

% Function for ODE system
function xd = OneDOFMSD_ex(t, x)
    k = 9.8;  % Spring constant (N/m^2)
    m = 1;    % Mass (kg)
    c = 0.5;  % Damping coefficient (N s/m)
    g = 9.8;  % Gravitational acceleration (m/s^2)

    % ODE system: dx1/dt = velocity, dx2/dt = acceleration due to forces
    xd(1) = x(2);  % Position (y) is the velocity
    xd(2) = (-k/m) * (x(1)) - (c/m) * (x(2)) - g;  % Velocity (v) is the acceleration due to spring, damper, and gravity
    xd = xd';  % Return as a column vector
end


%% Question 4: final code
% Initial condition
y0 = [0.5, 0];  
% Time range
time = 0:0.01:10;  

% Define the differential equation function
function diff_eq = ode_func(t, x)
    diff_eq = zeros(2,1);
    % The first equation is y' = v (velocity)
    diff_eq(1) = x(2);
    % The second equation: v' = acceleration
    diff_eq(2) = (+9.8 - 9.8 * (x(1)^2) - 0.5 * x(2));  
end

% Solve the system using ode45
[t, y] = ode45(@ode_func, time, y0);

% Plotting the position (y) versus time
figure;
plot(t, y(:,1));
xlabel('Time (s)');
ylabel('Position (y)');
title('Position vs Time for a 1D Damped Spring System');
grid on;
