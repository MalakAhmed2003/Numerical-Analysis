%%
clc
clear all
%% Question 1
% Part A
function [x_mat, temp_coefficient_mat, temp_sol_mat] = Gaussian_elimination_solver(coefficient_mat, sol_mat, n)
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

% Part B
coefficient_mat = [1 4 9 16;
4 9 16 25;
9 16 25 36;
16 25 36 49.0001];
sol_mat_b = [4; 400; 40000; 400000000];
n = length(sol_mat_b);
[x_matrix, new_coefficient_mat, temp_sol_mat] = Gaussian_elimination_solver(coefficient_mat, sol_mat_b, n);
disp('Updated Coefficient Matrix:');
disp(new_coefficient_mat);
disp('Updated Solution Vector:');
disp(temp_sol_mat);
disp('Solution Vector x:');
disp(x_matrix);

% verifying using the backslash method
solution_backslash = coefficient_mat\sol_mat_b

% Comparing the answers of my solver againist the backslash answer

x_relative_error = [];
for v = 1:1:n
    x_relative_error(v) = (x_matrix(v)- solution_backslash(v))/x_matrix(v);
end
disp('Error in each value of the calculated x_matrix:');
disp(x_relative_error)

%% Question 2
% Part a

function [L, U] = L_U_decomposer(coefficient_mat, n)
    if (det(coefficient_mat) == 0)
        disp('The system has no solution.');
        L = [];
        U = [];
        return;
    else
        U = coefficient_mat;
    L = eye(n);        
        % Gaussian elimination part
        for j = 1:n
            for i = (j+1):n
                % Calculating the elimination factor
                factor = U(i, j) / U(j, j);
                
                % Updatating the rows of the U matrix after
                % operating on them
                U(i, j:n) = U(i, j:n) - factor * U(j, j:n);
                
                % Updating the L matrix as well
                L(i, j) = factor;
            end
        end
       
    end
end

A =[
8 6 1;
10 6 14;
11 15 14];
b =[19; 18; 7];
[L, U] = L_U_decomposer(A, 3);
disp("U using the designed function")
disp(U)
disp("L using the designed function")
disp(L)
% Verifying that L*U = A
L*U

%% Question 2
% Part b
% LU_built_in checker

A =[
8 6 1;
10 6 14;
11 15 14];
[l,u, p] = lu(A);
disp("U using the built-in function")
disp(u)
disp("L using the built-in function")
disp(l)
% L*U must now equal A*P which is A multiplied permutation matrix that
% tracks row swaps done for pivoting
l*u
p*A
%% Question 2
% Part c
% Solving the values of X using the L_U solved matrices
A =[
8 6 1;
10 6 14;
11 15 14];
b =[19; 18; 7];
d = L\b;
x = U\d
% checking the value of x using the backslash operator
x_check = A\b
%% Question 3

M =[1 4 9 16 25;
    4 9 16 25 36;
    9 16 25 36 49;
    16 25 36 49 64;
    25 36 49 64 81]
M_inverse = inv(M);

% Part a
norm_M_l2 = norm(M, 2)
norm_M_inverse_l2 = norm(M_inverse, 2)
cond_num_M_l2 = norm_M_l2 * norm_M_inverse_l2


% Part b
norm_M_row_sum = 0;
norm_M_inverse_row_sum=0;
added_num = 0;
% Calculating the row sum norm
for i=1:1:5
    for j = 1:1:5
        if M(i,j)> added_num
             added_num= M(i,j);
        end
       
        
    end
     norm_M_row_sum = norm_M_row_sum + added_num;
     added_num = 0;
     
end

for l=1:1:5
    for o = 1:1:5
        if M_inverse(l,o)> added_num
             added_num= M_inverse(l,o);
        end
       
        
    end
     norm_M_inverse_row_sum = norm_M_inverse_row_sum + added_num;
     added_num = 0;
     
end
norm_M_row_sum
norm_M_inverse_row_sum

cond_num_M_row_sum = norm_M_row_sum * norm_M_inverse_row_sum


% Part c
lost_digits = round(log10(cond_num_M_l2))


% Part d
% Method 1: check whether A*A^-1 == I or not
M_M_inverse = M * M_inverse
 if M_M_inverse == eye(5)
     disp("The Matrix M is well-conditioned")
 else
     disp("The Matrix M is ill-conditioned")
 end

 % Method 2: check whether (A^-1)^-1 == A or not
inverse_M_inverse = inv(M_inverse)
 if inverse_M_inverse == M
     disp("The Matrix M is well-conditioned")
 else
     disp("The Matrix M is ill-conditioned")
 end
 