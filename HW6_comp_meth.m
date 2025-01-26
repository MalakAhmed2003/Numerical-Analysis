%%
clc 
clear all

%% Part a: 2-D regression using normal equations
data = load('C:\Users\DELL\Documents\MATLAB\HW6_PowerPlant_Data.mat');

disp(fieldnames(data));
humidity = data.humidity;
temp = data.temp;
power_output = data.power_output;
n = length(power_output);
average_power_output = mean(power_output); % Calculate average power output

% Initializing the coefficient matrix that we will use for the calculation
% of a_0, a_1, a_2 and then filling it according to the normal equations
A = ones(3, 3);
b = ones(3, 1);
A(1, 1) = n;
A(1, 2) = sum(humidity);
A(1, 3) = sum(temp);
A(2, 1) = sum(humidity);
A(2, 2) = sum(humidity .* humidity);
A(2, 3) = sum(humidity .* temp);
A(3, 1) = sum(temp);
A(3, 2) = sum(humidity .* temp);
A(3, 3) = sum(temp .* temp);
b(1, 1) = sum(power_output);
b(2, 1) = sum(humidity .* power_output);
b(3, 1) = sum(temp .* power_output);

% Function for Gaussian elimination solver
function [x_mat] = Gaussian_elimination_solver(coefficient_mat, sol_mat)
    if (det(coefficient_mat) == 0)
        disp('The system has no solution.');
        x_mat = [];
        return;
    else 
        n = length(sol_mat);
        x_mat = zeros(n, 1);
        temp_coefficient_mat = coefficient_mat;
        temp_sol_mat = sol_mat;

        % Gaussian elimination part
        for j = 1:n
            for i = (j+1):n
                % Calculating the elimination factor
                factor = temp_coefficient_mat(i, j) / temp_coefficient_mat(j, j);
                
                % Updating the rows of the coefficient matrix
                temp_coefficient_mat(i, j:n) = temp_coefficient_mat(i, j:n) - factor * temp_coefficient_mat(j, j:n);
                
                % Updating the solution matrix
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

% Call the Gaussian elimination solver
a = Gaussian_elimination_solver(A, b);


% Plotting of the data points and the meshgrid representing the surface
% that best fits the data points
figure1 = figure;
axes1 = axes('Parent', figure1);
hold(axes1, 'on');

% 3D plot of data using the independent variables (humidity and
% temperature) and the dependent variable power_output
plot3(humidity, temp, power_output, 'MarkerFaceColor', [1 0 1], 'MarkerSize', 20, 'Marker', '.', 'LineStyle', 'none', 'Color', [1 0 1]);

% Setting the labels of the axes
xlabel('Humidity');
ylabel('Temp');
zlabel('Power Output');

% Creating the surface for the regression plane that best fits the data
% points according to the calculated coefficients using the normal
% equations method
[humidity_mesh, temp_mesh] = meshgrid(linspace(min(humidity), max(humidity), 10), linspace(min(temp), max(temp), 10));
z = a(1) + a(2) * humidity_mesh + a(3) * temp_mesh;
surf(humidity_mesh, temp_mesh, z, 'Parent', axes1, 'MarkerEdgeColor', 'none', 'FaceColor', 'none', 'EdgeColor', [0.39215686917305 0.474509805440903 0.635294139385223]);
view(axes1, [-9.97999999999995 41.52]);


% r², s_r, s_y_x calculations
residuals = power_output - (a(1) + a(2) * humidity + a(3) * temp);
s_r = sum(residuals.^2) % Calculate residual sum of squares

s_y_x = sqrt(s_r / (n - 2)) % Calculation of the standard error of the estimate

s_t = sum((power_output - average_power_output).^2);
r_2 = (s_t - s_r) / s_t 

%% Part b: 2_D regression using matlab method

% Loading the data from the file to arrays representing the independent
% variables and the output
data = load('C:\Users\DELL\Documents\MATLAB\HW6_PowerPlant_Data.mat');
humidity = data.humidity;
temp = data.temp;
power_output = data.power_output;
n = length(power_output);
average_power_output = sum(power_output) / n;
first_column = ones(n, 1);
A_ML = [first_column, humidity, temp];
a_ML = A_ML\power_output;

% Plotting of the data points and the meshgrid representing the surface
% that best fits the data points
figure2 = figure;
axes1 = axes('Parent', figure2);
hold(axes1, 'on');

% 3D plot of data using the independent variables (humidity and
% temperature) and the dependent variable power_output
plot3(humidity, temp, power_output, 'MarkerFaceColor', [1 1 0], 'MarkerSize', 20, 'Marker', '.', 'LineStyle', 'none', 'Color', [1 0 0]);

% Setting the labels of the axes
xlabel('Humidity');
ylabel('Temp');
zlabel('Power Output');

% Creating the surface for the regression plane that best fits the data
% points according to the calculated coefficients using the MATLAB direct
% method
[humidity_mesh, temp_mesh] = meshgrid(linspace(min(humidity), max(humidity), 10), linspace(min(temp), max(temp), 10));
z = a_ML(1) + a_ML(2) * humidity_mesh + a_ML(3) * temp_mesh;
surf(humidity_mesh, temp_mesh, z, 'Parent', axes1, 'MarkerEdgeColor', 'none', 'FaceColor', 'none', 'EdgeColor', [0.39215686917305 0.474509805440903 0.635294139385223]);
view(axes1, [-9.97999999999995 41.52]);

% r², s_r, s_y_x calculations
residuals_ML = power_output - (a_ML(1) + a_ML(2) *humidity + a_ML(3) * temp);
s_r_ML = sum(residuals_ML.^2)

s_y_x_ML = sqrt(s_r_ML / (n - 2)) % Calculation of the standard error of the estimate

s_t_ML = sum((power_output - average_power_output).^2);
r_2_ML = (s_t_ML - s_r_ML) / s_t_ML


%% Part c: 2_D regression using the “inverse of the matrix” method
% Loading the data from the file to arrays representing the independent
% variables and the output
data = load('C:\Users\DELL\Documents\MATLAB\HW6_PowerPlant_Data.mat');

humidity = data.humidity;
temp = data.temp;
power_output = data.power_output;

n = length(power_output);
I = ones(n, 1);
average_power_output = sum(power_output) / n;
first_column = ones(n, 1);
A_inverse_method = cat(2,first_column, humidity, temp); % Concatinating the variables to form the A matrix fr the system as Aa = b

% Calculate coefficients using the inverse method
a_inverse_method = inv((A_inverse_method') * A_inverse_method) * (A_inverse_method') * power_output;

% Plotting of the data points and the meshgrid representing the surface
% that best fits the data points
figure3 = figure;
axes1 = axes('Parent', figure3);
hold(axes1, 'on');

% 3D plot of data  using the independent variables (humidity and
% temperature) and the dependant variable power_output
plot3(humidity, temp, power_output, 'MarkerFaceColor', [1 1 0], 'MarkerSize', 20, 'Marker', '.', 'LineStyle', 'none', 'Color', [0 0 1]);

% Setting the labels of the axes
xlabel('Humidity');
ylabel('Temp');
zlabel('Power Output');


% Creating the surface for the regression plane that best fits the data
% points according to the calculated coefficients using the inverse method
[humidity_mesh, temp_mesh] = meshgrid(linspace(min(humidity), max(humidity), 10), linspace(min(temp), max(temp), 10));
z = a_inverse_method(1) + a_inverse_method(2) * humidity_mesh + a_inverse_method(3) * temp_mesh;
surf(humidity_mesh, temp_mesh, z, 'Parent', axes1, 'MarkerEdgeColor', 'none', 'FaceColor', 'none', 'EdgeColor', [0.39215686917305 0.474509805440903 0.635294139385223]);
view(axes1,[-9.97999999999995 41.52]);

% r², s_r, s_y_x calculations
residuals_inverse_method = power_output - (a_inverse_method(1) + a_inverse_method(2) * humidity + a_inverse_method(3) * temp);
s_r_inverse_method = sum(residuals_inverse_method.^2)

s_y_x_inverse_method = sqrt(s_r_inverse_method / (n - 2)) % Calculation of the standard error of the estimate

s_t_inverse_method = sum((power_output - average_power_output).^2);
r_2_inverse_method = (s_t_inverse_method - s_r_inverse_method) / s_t_inverse_method

