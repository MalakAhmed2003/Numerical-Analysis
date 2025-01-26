%%
clc
clear all
%% Question 1
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
