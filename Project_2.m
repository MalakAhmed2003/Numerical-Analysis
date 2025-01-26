%%
clc
clear all
%%Chapter 4: part1
%{ hh}%
a = -3;
h = -1.5;
error = zeros(1,4);
f = @(x) -1*25*(x^3) + 8*(x^2) + 7*x + 20;
f_1 = @(x) -1*75*(x^2) + 16*x + 7; %first derivative
f_2 = @(x) -1*150*(x) + 16;        %second derivative
f_3 = @(x) -150;                   %third derivative
f_order = zeros(1,4);              %function up to nth order
% calculating the value of the function using taylor expansion, order 0
f_order(1) = f(a); 
% calculating the value of the function using taylor expansion, order 1
f_order(2) = f(a) + (f_1(a)*h);  
% calculating the value of the function using taylor expansion, order 2
f_order(3) = f(a) + (f_1(a)*h) + (f_2(a)*(h^2)/factorial(2));
% calculating the value of the function using taylor expansion, order 3
f_order(4) = f(a) + (f_1(a)*h) + (f_2(a)*(h^2)/factorial(2)) + (f_3(a)*(h^3)/factorial(3));  
for i=1:1:4
  error(i) = ((f(-4.5)-f_order(i))/f(-4.5))*100;
end
for i=1:1:4
   fprintf('For order = %f', i-1)
   fprintf(' The value of f(-4.5) = %f', f_order(i))
   fprintf(' The error = %f', error(i))
   fprintf('\n ')
end

% Data display in a table
data = {
    "0", f_order(1), error(1);
    "1", f_order(2), error(2);
    "2", f_order(3), error(3);
    "3", f_order(4), error(4);
};
T = cell2table(data, 'VariableNames', {'Order', 'f(-4.5)', 'Relative True Error'});
disp(T);
%% Chapter 4: Part 2
f = @(x) -1*25*(x^3) + 8*(x^2) + 7*x + 20;
f_1 = @(x) -1*75*(x^2) + 16*x + 7; %first derivative
h = .5;
x_i = -4.5;
x_i_plusone = x_i + h;
x_i_minusone = x_i - h;
% using forward difference
f_forward_derv = (f(x_i_plusone)-f(x_i))/h

error_forward = ((f_1(x_i)-f_forward_derv)/f_1(x_i))*100

% using forward difference
f_backward_derv = (f(x_i)-f(x_i_minusone))/h

error_backward = ((f_1(x_i)-f_backward_derv)/f_1(x_i))*100

% using forward difference
f_centr_derv = (f(x_i_plusone)-f(x_i_minusone))/(2*h)

error_centr = ((f_1(x_i)-f_centr_derv)/f_1(x_i))*100

% Data display in a table
data = {
    "Forward Difference", f_forward_derv, error_forward;
    "Backward Difference", f_backward_derv, error_backward;
    "Central Difference", f_centr_derv, error_centr;
    };
T = cell2table(data, 'VariableNames', {'Method Used', 'first derivative at -4.5', 'Relative True Error'});
disp(T);
%% Chapter 4: part 3
delta_c = 3;
% defining the function of velocity in terms of m, g, t and c
syms t c m g
v = (g*m/c)*(1 -exp(-1*c*t/m));
% differentiating the function with respect to c
dv_by_dc = diff(v, c);
% substituting in the derivative function by the values of m, g, t and c
abs_dv_by_dc = abs(double(subs(dv_by_dc, {g, m, c, t}, {9.8, 50, 9.5, 4})));
% calcultating the error in v
delta_v = delta_c * abs_dv_by_dc;
fprintf('The error of v = %f', abs_dv_by_dc)

%% Chapter 4: Part 4

x_a = 1.0001;
x_b = 9;
x_c = 200;
x_d = 0.01;

syms x
% Defining the function in part a, its derivative and calculating its
% condition number
f_a = @(x) sqrt(abs(x-1))+1;

dirv_f_a = diff(f_a,x);
dirv_x_a = subs(dirv_f_a, {x}, {x_a});
cond_num_a = double(dirv_x_a*x_a/f_a(x_a))

% Defining the function in part b, its derivative and calculating its
% condition number
f_b = @(x) exp(-1*x);

dirv_f_b = diff(f_b,x);
dirv_x_b = subs(dirv_f_b, {x}, {x_b})
cond_num_b = double(dirv_x_b*x_b/f_b(x_b))

% Defining the function in part c, its derivative and calculating its
% condition number
f_c = @(x) (sqrt((x^2)+1))-x;

dirv_f_c = diff(f_c, x);
dirv_x_c = subs(dirv_f_c, {x}, {x_c})
cond_num_c = double(dirv_x_c*x_c/f_c(x_c))

% Defining the function in part d, its derivative and calculating its
% condition number
f_d = @(x) (exp(x)-1)/x;

dirv_f_d = diff(f_d, x);
dirv_x_d = subs(dirv_f_d, {x}, {x_d})
cond_num_d = double(dirv_x_d*x_d/f_d(x_d))

% Data display in a table
data = {
    "(a)", cond_num_a, "ill_conditioned >> 1";
    "(b)", cond_num_b, "moderately-well-conditioned";
    "(c)", cond_num_c, "well-conditioned = 1";
    "(d)", cond_num_d, "very-well-conditioned < 1";
};
T = cell2table(data, 'VariableNames', {'Function Number', 'Condition Number', 'Interpretation'});
disp(T);

%% Chapter 5: part a
% Define the range and function
x = -100:1:100; 
f = @(x) 10 + 45.4*x - 10*(x.^2) + 45.4*(x.^3) - 9.7*(x.^4) + 4*(x.^5); 

% Plot the function
fplot(f, [-100, 100]); 
xlabel('x'); 
ylabel('f(x)'); 
title('Plot of the function f(x)'); 
grid on; 

% Finding the roots 
roots = fzero(f, [-100, 100])
% Marking the intersection points on the graph
hold on;
plot(roots, f(roots), 'ro', 'MarkerSize', 8, 'DisplayName', 'Intersection Points'); % Red circles for intersections
legend show; 
hold off; 

%% Chapter 5: part b
% Define the function
f = @(x) 10 + 45.4*x - 10*(x.^2) + 45.4*(x.^3) - 9.7*(x.^4) + 4*(x.^5);

% Define the bisection function
function [root, iter, results] = bisection_method(f, x_l, x_u, criterion_error, max_iter)
    relative_error = 1000;
    iter = 1;
    c = (x_l + x_u) / 2; % Initial guess
    results = []; % Initialize results array

    while (relative_error > criterion_error && iter < max_iter)
        c_old = c; % Store the old value of c
        c = (x_l + x_u) / 2; % Calculate c

        % Calculate relative error after the first iteration
        if iter > 1
            relative_error = abs((c - c_old) / c) * 100;
        end

        results = [results; iter, x_l, x_u, c]; % storing the values of the current iteration in the results array

        % Update the bounds of the bracketing interval
        if (f(c) * f(x_l) < 0)
            x_u = c; 
        else
            x_l = c; 
        end
                
        root = c; 
        iter = iter + 1; % Increment iteration counter
    end
end

% Call the function
format long
[root, iter, results] = bisection_method(f, -0.5, 0, 10, 1000);

% Create a table to display the results
T = array2table(results, 'VariableNames', {'Iteration', 'x_l', 'x_u', 'Root'});

% Display the table
disp(T);

%% Chapter 5: part c

% Define the function
f = @(x) 10 + 45.4*x - 10*(x.^2) + 45.4*(x.^3) - 9.7*(x.^4) + 4*(x.^5);

% Define the false position function
function [root, iter, results] = false_position(f, x_l, x_u, criterion_error, max_iter)
    relative_error = 1000;
    iter = 1;
     c = ((-1 * f(x_l) * (x_u - x_l)) / (f(x_u) - f(x_l))) + x_l; % Initial guess
     results = []; % Initialize results array

    while (relative_error > criterion_error && iter < max_iter)
        c_old = c; % Store the old value of c
        % Calculate c using the false position formula
        c = ((-1 * f(x_l) * (x_u - x_l)) / (f(x_u) - f(x_l))) + x_l

        results = [results; iter, x_l, x_u, c]; % storing the values of the current iteration in the results array

        % Update the bounds of the bracketing interval
        if (f(c) * f(x_l) < 0)
            x_u = c 
        else
            x_l = c 
        end

        % Calculate relative error after  the first iteration
        if iter > 1
            relative_error = abs((c - c_old) / c) * 100
        end
        
        root = c; 
        iter = iter + 1 % Increment iteration counter
    end
end


% Call the function
format long
[root, iter, results] = false_position(f, -0.5, 0, .1, 1000)

% Create a table to display the results
T = array2table(results, 'VariableNames', {'Iteration', 'x_l', 'x_u', 'Root'});

% Display the table
disp(T);
