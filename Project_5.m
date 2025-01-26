%% Data
%%  PreferVegemite CityPopulation
             194300 1815000
             192700 1805000
             193500 1795000
             179900 1793000
             179100 1765000
             180700 1755000
             183100 1740000
             177500 1740000
             176700 1735000
             176220 1720000
             175100 1725000
             165500 1715000
             163900 1710000
             163900 1710000
             163100 1705000
             161500 1645000
             151100 1655000
             117996 1650000
             149500 1645000
             152700 1635000
             154300 1630000
             138300 1620000
             128540 1605000
             126188 1610000
             128380 1605000
             153980 1615000
             155100 1625000
clc
clear all
%% Question 1

% Problem: Vector and Matrix Norm Calculations

% For the following vectors and matrices (Note: You may use MATLAB's built-in "norm" functions for Question 1):

% V1 = [-1.3 1.4 3.2];
% V2 = [2.8 2.4 3.1];

% M1 = [2.3 2.2 7.3;
%       3.9 7.9 -3.3;
%       7.7 3.5 -0.5];

% M2 = [-1.6 2.8 9.2;
%       2.5 7.9 -1.3;
%       6.4 3.3 -0.5];

% Part (a):
% Calculate the difference between the two vectors, V1 and V2, and quantify their magnitude using the following norm measures.

% i. Calculate the L1 norm of the vector difference.
% ii. Calculate the L2 norm of the vector difference.
% iii. Calculate the Infinity norm of the vector difference.

% Part (b):
% Calculate the difference between the two matrices, M1 and M2, and quantify their magnitude using the following norm measures.

% i. Calculate the L1 norm of the matrix difference.
% ii. Calculate the L2 norm of the matrix difference.
% iii. Calculate the Frobenius norm of the matrix difference.
% iv. Calculate the Infinity norm of the matrix difference.


% Question 1
% Part a
V_1 = [-1.3, 1.4, 3.2];
V_2 = [2.8, 2.4, 3.1];
V = V_1 - V_2; % Finding the difference between the two vectors
L_1_norm = norm(V, 1);
L_2_norm = norm(V, 2);
L_inf_norm = norm(V, Inf);
data = {
    "(L1)", L_1_norm, "L1 norm returns the sum of the absolute values of the elements, soit gives the greatest value compared to other norms.";
    "(L2)", L_2_norm, "L2 norm returns the square root of the sum of the square of elements in the vector";
    "(Infinity Norm)", L_inf_norm, "The infinity norm returns the largest absolute value of the elements in the vector, which returns the smallest value.";
   
};
T = cell2table(data, 'VariableNames', {'Norm Type', 'Norm Value', 'Comments'});
disp(T);


%% Question 1
% Part b

M_1 =[2.3 2.2 7.3;
3.9 7.9 -3.3;
7.7 3.5 -0.5];
M_2 =[-1.6 2.8 9.2;
2.5 7.9 -1.3;
6.4 3.3 -0.5];
M = M_1 - M_2

L_1_norm = norm(M, 1);
L_2_norm = norm(M, 2);
L_fro_norm = norm(M, 'fro');
L_inf_norm = norm(M, Inf);

data = {
    "L1", L_1_norm, "This L1 returns the maximum of the column sum in the matrix.";
    "L2", L_2_norm, "This L2 depends on calculating the square root of the maximum eigenvalues of the matrix.";
    "Frobenius Norm", L_fro_norm, "This norm returns the square root of the sum of the squares of all the elements of the matrix.";;
    "Infinity Norm", L_inf_norm, "This Infinity norm returns the maximum of the row sum of absolute values.";
   
};
T = cell2table(data, 'VariableNames', {'Norm Type', 'Norm Value', 'Comments'});
disp(T);

%% Question 2

% Problem: Gauss-Seidel Method for Solving NxN System of Equations

% Part (a):
% Write a MATLAB program to solve an NxN system of equations of the form A * x = b, 
% implementing the Gauss-Seidel algorithm. 
% Use a stopping criterion such that the worst-case (the vector component with the largest approximate error) 
% answer is correct to 3 significant digits. 
% Ensure that you implement a maximum iteration count in your loop, check for the Gauss-Seidel convergence criteria.
% Part (b):
% Use the code you wrote to compute the solution, and then check your solution by computing 
% the solution using MATLAB’s backslash operator. 
% Show the error between the true solution and the solution obtained from the Gauss-Seidel method.

% Given matrix A:
% A = [20 4 2 3;
%      1 34 8 6;
%      3 4 18 2;
%      9 8 5 69];

% Given vector b:
% b = [5;
%      2;
%     -7;
%     12];



function [x] = Gauss_seidel_solver(A, b)
    n = length(b);
    x = zeros(n, 1); % Initializing the solutions as zeros
    max_error = 1000;  % Initialize max_error to a large number for the first iteration at least to run
    stopping_error = 0.5 * (10^(2-3));  % Stopping criterion
    max_iterations = 1000;  % Maximum number of iterations
    iter = 0;
    
    diagonal_dominance = true;
    for k=1:1:n
        summation = 0;
        for q=1:1:n
            if q ~= k
                summation = summation + abs(A(k,q));
            end
        end

        if A(k,k) > summation % Checking the diagonal dominance f the matrix
        else
            diagonal_dominance = false;
        end
    end
    if diagonal_dominance == true
        disp('Gauss-siedel method is guaranteed to converge, since the matrix is diagonally dominant.');
        while iter < max_iterations && max_error > stopping_error
        x_old = x;  % Store the last estimate in the x_old array
        for i = 1:n
            sum = b(i);
            for j = 1:n
                if j ~= i
                    sum = sum - A(i, j) * x(j); % Calculating the new estimate of x as defined by the Gauss siedel method
                end
            end
            x(i) = sum / A(i, i); 
        end
        
        % Calculating the maximum error
        max_error = max(abs((x - x_old) ./ (x)) * 100);  
        iter = iter + 1;  % Increment iteration counter
    end
    else
        disp('Gauss-siedel method is not guaranteed to converge, since the matrix is not diagonally dominant.');
    end
    
end

% Part b
A =[
20 4 2 3;
1 34 8 6;
3 4 18 2;
9 8 5 69]
b =[
5;
2;
-7;
12]

x_using_Gauss_solver = Gauss_seidel_solver(A, b);
disp("The solutions using Gauss-siedel solver are:")
disp(x_using_Gauss_solver);
% Checking the output solution 
x_backslash = A\b;
disp("The solutions using Backslash are:")
disp(x_backslash);
Error_between_the_backslash_and_our_solver = (x_backslash - x_using_Gauss_solver) % showing the error between backslash method and Gauss-siedel, which is actually almost zero


%% Question 3

% Problem: Analyzing Population and Vegemite Preference Data

% You have data on the population of several cities, 
% and the number of individuals from each city who prefer Vegemite over other spreadable food products.

% Part (a):
% Write a MATLAB program to calculate the equation of the best-fitting straight line 
% (1D linear regression) for the given X and Y data using the "Normal Equations" method.

% Part (b):
% Generate a plot that displays both the original data points and the computed regression line on the same graph, 
% allowing you to visually assess how well the regression line fits the data.

% Part (c):
% Calculate and print the standard error of the estimate (Sy/x), and provide the values for the regression coefficients a0 and a1.



% Loading the data
data_loaded = importdata('HW5_P3data.mat');
n = length(data_loaded.x);
I = ones(27,1);

y_sum = sum(data_loaded.y);
x_sum = sum(data_loaded.x);

x_2_sum = sum((data_loaded.x).*(data_loaded.x));
x_y_sum = sum((data_loaded.x).*(data_loaded.y));
a_1 = (n*x_y_sum - x_sum*y_sum)/(n*x_2_sum - ((x_sum)^2))
a_0 = (y_sum - a_1*x_sum)/n

% Part b
% Plotting the scattered data points and the line we found to see 
hold on;
plot(data_loaded.x, data_loaded.y, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
 
f = a_0 +a_1 * data_loaded.x; % Definition of the straight line that best fit the data after finding a_0 and a_1
plot(data_loaded.x, f, 'b')
xlabel('City Population');
ylabel('Prefer Vegemite');

% Part c

s_r = sum((data_loaded.y -a_1*data_loaded.x - a_0*I ).*(data_loaded.y -a_1*data_loaded.x - a_0*I )) % sum of quares calculation

s_y_x = sqrt(s_r/(n-2)) % calculating the standard error of the estimate
