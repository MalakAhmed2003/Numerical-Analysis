%% 
clear all
clc
%% Question 1

% Find the Highest Real Root of f(x) = x^3 - 13.1x^2 + 47x - 34.8:
% 1. Graphical Method: Use MATLAB to plot the function and visually identify the highest real root. 
% 2. Fixed-Point Iteration Method: Apply the Fixed-Point Iteration method until the approximation error is less than 0.1%. Start with an initial value of x_0 = 8. 
% 3. Newton-Raphson Method: Use the Newton-Raphson method to find the root, with a maximum of 6 iterations and an initial guess of x_0 = 8. 
% 4. Secant Method: Apply the Secant method to find the root, with a maximum of 6 iterations. Use the initial guesses x−1=20x_{-1} = 20x−1​=20 and x_0 = 8.

%%part 1
% Define the function
f = @(x) (x.^3) - 13.1*(x.^2) + 47*x - 34.8;

% Calculate the roots
roots_f = roots([1, -13.1, 47, -34.8]);

% Create a range of x values for plotting
x = linspace(-1, 10, 400); 
y = f(x);
% Plot the function
figure;
plot(x, y, 'b-', 'LineWidth', 2); % Plot the function in blue
ylim([min(f(x)), max(f(x))])
hold on;

% Plot the roots
plot(roots_f, f(roots_f), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

% Add grid and labels
grid on;
xlabel('x');
ylabel('f(x)');
title('Plot of the function f(x) and its roots');
legend('f(x)', 'Roots', 'Location', 'Best');

hold off;
disp(max(roots_f))
%% Question 1: part 2

% Define the function
g = @(x)  ((x.^3) - 13.1*(x.^2) - 34.8)/-47;
x_i = 8;
relative_err = 100;
max_iter = 1000;
iter = 0;
results = [];
while (relative_err > 0.1 && iter < max_iter)
    x_old = x_i;
    x_i = g(x_i) %Doing the fixed point update
    
    % Calculate the relative error with respect to the previous iteration
    relative_err = abs((x_i - x_old) / x_i) * 100
    
    iter = iter + 1
    results = [results; iter, x_i, relative_err]
end
T = array2table(results, 'VariableNames',{'Iteration number(i)', 'x_i', 'Relative error'});
disp(T)
%% Question 1: part 3

% Define the function
g = @(x)  ((x.^3) - 13.1*(x.^2) + 47*x - 34.8);

% Define the derivative function using symbolic differentiation
syms x_sym
g_sym = ((x_sym.^3) - 13.1*(x_sym.^2) + 47*x_sym - 34.8);
dirv_g = matlabFunction(diff(g_sym, x_sym));  % Convert to function handle

x_i = 8;  % Initial guess
max_iter = 6;
iter = 0;
results = [];

while (iter < max_iter)
    x_old = x_i;
    x_i = x_i - (g(x_i) / dirv_g(x_i));  % Newton's method update
    
   
    iter = iter + 1;
    
    % Append results
    results = [results; iter, x_i];
end

% Create a table to display results
T = array2table(results, 'VariableNames', {'Iteration number(i)', 'x_i'});
disp(T);

%% Question 1: part 4

% Define the function
g = @(x)  ((x.^3) - 13.1*(x.^2) + 47*x - 34.8);

% Define the derivative function using symbolic differentiation
syms x_sym
g_sym = ((x_sym.^3) - 13.1*(x_sym.^2) + 47*x_sym - 34.8);


x_0 = 20
x_i = 8;  % Initial guess
max_iter = 6;
iter = 0;
results = [];

while (iter < max_iter)
    x_old = x_i;
    x_i = x_i - (g(x_i) / ((g(x_i)-g(x_0))/(x_i - x_0)));  % Secant method update
    x_0 = x_old; 
       
    iter = iter + 1;
    
    % Append results
    results = [results; iter, x_i];
end

% Create a table to display results
T = array2table(results, 'VariableNames', {'Iteration number(i)', 'x_i'});
disp(T);

%% Question 2

% 1. **Gaussian Elimination Solver**: Write a MATLAB function to solve a system of linear equations using the Gaussian Elimination method. Your solver should be able to handle both a 3x3 matrix and a general NxN system of equations. 
% **Note**: In Homework 3, you have the option to implement the solver for just a 3x3 matrix, but in Homework 4, you'll be required to extend it to an NxN system. 
% Writing the solver for the NxN system now will help you get ahead for the next assignment.

% 2. **Apply the Solver**: Use the Gaussian Elimination code you developed in part 1 to solve the system of equations:
% |6 10 16  | | x_1 |    | 6   |
% |-18 10 24| | x_2 | =  | -30 |
% |24 50 60 | | x_3 |    | 38  |
% The system is in the form Ax = b, where you need to find x. 

%Part B :1
% 3*3 system of linear equations solver(Gaussian Elemination method)

function [x_mat] = Guassian_elemination_solver(coefficient_mat, sol_mat)
if (det(coefficient_mat) == 0)
    disp('The system has no solution.')
else 
    x_mat = zeros(3,1);
    temp_coefficient_mat = coefficient_mat;
    temp_sol_mat = sol_mat;
    temp_sol_mat(2) = sol_mat(2) - coefficient_mat(2,1) * sol_mat(1)/ coefficient_mat(1, 1);
    temp_coefficient_mat(2, [1:3]) = coefficient_mat(2, [1:3]) - coefficient_mat(2,1) * coefficient_mat(1, [1:3])/ coefficient_mat(1, 1);
    temp_sol_mat(3) = sol_mat(3) - coefficient_mat(3,1) * sol_mat(1)/ coefficient_mat(1, 1);
    temp_coefficient_mat(3, [1:3]) = coefficient_mat(3, [1:3]) - coefficient_mat(3,1) * coefficient_mat(1, [1:3])/ coefficient_mat(1, 1);
    temp_sol_mat(3) = temp_sol_mat(3) - temp_coefficient_mat(3,2) * temp_sol_mat(2)/ temp_coefficient_mat(2, 2);
    temp_coefficient_mat(3, [1:3]) = temp_coefficient_mat(3, [1:3]) - temp_coefficient_mat(3,2) * temp_coefficient_mat(2, [1:3])/ temp_coefficient_mat(2, 2);
    x_mat(3) = temp_sol_mat(3)/ temp_coefficient_mat(3,3);
    x_mat(2) = (temp_sol_mat(2)- (x_mat(3)*temp_coefficient_mat(2,3)))/temp_coefficient_mat(2,2);
    x_mat(1) = (temp_sol_mat(1)- (x_mat(3)*temp_coefficient_mat(1,3)) - (x_mat(2)*temp_coefficient_mat(1,2)))/temp_coefficient_mat(1,1);
end
end


% n*n system of linear equations solver(Gaussian Elemination method)
function [x_mat, temp_coefficient_mat, temp_sol_mat] = Gaussian_elimination_solver_n(coefficient_mat, sol_mat, n)
    if (det(coefficient_mat) == 0)
        disp('The system has no solution.');
        x_mat = [];
        temp_coefficient_mat = [];
        temp_sol_mat = [];
        return;
    else 
        x_mat = zeros(n, 1);
        temp_coefficient_mat = coefficient_mat;
        temp_sol_mat = sol_mat;
        
        % Gaussian elimination part
        for j = 1:n
            for i = (j+1):n
                % Calculating the elimination factor
                factor = temp_coefficient_mat(i, j) / temp_coefficient_mat(j, j);
                
                % Updatating the rows of the coefficient matrix after
                % operating on them
                temp_coefficient_mat(i, j:n) = temp_coefficient_mat(i, j:n) - factor * temp_coefficient_mat(j, j:n);
                
                % Updating the solution matrix as well
                temp_sol_mat(i) = temp_sol_mat(i) - factor * temp_sol_mat(j);
            end
        end
        
        % Back substitution part
        x_mat(n) = temp_sol_mat(n) / temp_coefficient_mat(n, n);
        for k = n-1:-1:1
            x_mat(k) = temp_sol_mat(k);
            for z = k+1:n
                x_mat(k) = x_mat(k) - (x_mat(z) * temp_coefficient_mat(k, z));
            end
            x_mat(k) = x_mat(k) / temp_coefficient_mat(k, k);
        end
    end
end

% Part B:2
%sol_matrix_b refers to the b matrix in the right hand side of the
%equation, coefficient_mat refers to the A matrix in the left hand side of
%the equation
coefficient_mat = [6 10 16 ; -18 10 24; 24 50 60];
sol_mat_b = [6;-30;38];
x_matrix = Guassian_elemination_solver(coefficient_mat, sol_mat_b);
disp('solution matrix =')
disp(solution_matrix)
% comparing the values we got quantitatively by trying bcksubstitution in
% the original equation and comparing results
b_matrix_check = coefficient_mat*x_matrix


%% Question 3

% A Hilbert Matrix is a type of square matrix that is known for being ill-conditioned. 
% Each element Hij in a Hilbert Matrix is defined as:
% Hij = 1 / (i + j - 1)

% a) Generate a 10x10 Hilbert matrix and a 5x5 Hilbert matrix. 
% You can use a for-loop or a built-in MATLAB function to create them. 
% Name the matrices H10 and H5, respectively, and display them in your solution. 

% b) Create the column vectors X10_actual and X5_actual:
% X10_actual = [1; 2; 3; ...; 10]
% X5_actual = [1; 2; ...; 5]
% You can easily generate these vectors using shorthand notation [1:n]’.

% c) Using the Hilbert matrices H10 and H5, compute the column vectors:
% B10 = H10 * X10_actual
% B5 = H5 * X5_actual
% Then display the values of B10 and B5.

% d) Compute the condition numbers for both Hilbert matrices (H10 and H5) 
% and display the results.

% e) With the Hilbert matrices H10 and H5 and the column vectors B10 and B5, 
% solve for X using two methods: the backslash operator and the matrix inverse.
% Solve for X using the backslash operator (X_backslash) and also by using the inverse of the matrix (X_inverse).

% f) Calculate and display the errors for both methods and for both Hilbert matrices. 
% Use the following formulas for the errors:
% norm(X_actual - X_from_backslash) and norm(X_actual - X_from_inverse).
% Compare the results for both methods and matrices. 
% Determine which method is more accurate and explain why one matrix produces smaller errors.

% g) Using tic/toc, measure the time taken by your computer to solve the system 
% using each of the two methods and compare the timings. 
% NOTE: When measuring execution time, it’s recommended to use a for-loop with enough iterations 
% to make the execution last several tens of milliseconds. This helps ensure more accurate timing 
% by averaging out any small fluctuations caused by factors like disk access and pre-computation.


% Part a

% Define a function that creates the ill-conditioned Matrices
function [mat] = Hilbert_Matrix_creator(n)
    mat = zeros(n);
    for i=1:1:n
        for j= 1:1:n
            mat (i,j) = 1/(i + j - 1);

        end

    end
end
% Creating H10 and H5 matrices
H5 = Hilbert_Matrix_creator(5);
H10 = Hilbert_Matrix_creator(10);
% displaying the two matrices up to 2 decimal places
disp(num2str(H5, 2))
disp(num2str(H10, 2))
% Question 3: Part b
% Generate the column vectors X10actual and X5actual 
function [mat]= create_column_vector(n)
    mat = 1:1:n;
    mat= transpose(mat);
end
X10actual = create_column_vector(10)
X5actual = create_column_vector(5)

% Question 3: Part c  
B5 = H5 * X5actual
B10 = H10 * X10actual

% Question 3: Part d
%calculating the condition numbers for H5 and H10
 cond_H5 = cond(H5)
 cond_H10 = cond(H10)

% Question 3: Part e

% Solve for X from backslash operator
X_from_backslash_5 = H5 \ B5
X_from_backslash_10 = H10 \ B10

% Solve for X from inverse operator
X_from_inverse_5 = inv(H5) * B5
X_from_inverse_10 = inv(H10) * B10

% Question 3: Part f

% Error in X5actual calculation
err_backslash_5 = norm(X5actual - X_from_backslash_5)
err_inverse_5 = norm(X5actual - X_from_inverse_5)

% Error in X5actual calculation
err_backslash_10 = norm(X10actual - X_from_backslash_10)
err_inverse_10 = norm(X10actual - X_from_inverse_10)

% Question 3: Part g
% Calculating the time taken by each of the two methods
 num_iterations = 10;
 tic;
 for z=1:num_iterations
     X_from_backslash_5 = H5 \ B5;
 end
 time_backslash_5 = toc/num_iterations


 tic;
 for z=1:num_iterations
     X_from_inverse_5 = inv(H5) * B5;
 end
 time_inverse_5 = toc/num_iterations

 tic;
 for z=1:num_iterations
     X_from_backslash_10 = H10 \ B10;
 end
 time_backslash_10 = toc/num_iterations

 tic;
 for z=1:num_iterations
     X_from_inverse_10 = inv(H10) * B10;
 end
 time_inverse_10 = toc/num_iterations
