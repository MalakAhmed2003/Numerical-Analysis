clc
clear
%% Question 1
%{
Write a MATLAB program to calculate the roots of a quadratic equation of the form ax2+bx+c=0ax^2 + bx + c = 0ax2+bx+c=0. Implement this as a function, and use a subroutine (a sub-function) to compute the roots, which could be either real or complex.
First, design the algorithm by creating pseudocode. Then, implement the pseudocode in MATLAB. For testing, perform the following cases:
(a) a=2,b=7,c=4a = 2, b = 7, c = 4a=2,b=7,c=4
(b) a=0,b=−5,c=2.5a = 0, b = -5, c = 2.5a=0,b=−5,c=2.5
(c) a=2,b=1,c=8a = 2, b = 1, c = 8a=2,b=1,c=8
In your implementation, you are allowed to use basic MATLAB functions and your understanding of solving quadratic equations, but refrain from using MATLAB's built-in functions for solving equations. Additionally, ensure that you debug and document your code clearly.

}%
% Defining the values of the quadratic equation coefficients

a1 = 2;
b1 = 7;
c1 = 4;
a2 = 0;
b2 = -5;
c2 = 2.5;
a3 = 2;
b3 = 1;
c3 = 8;
% Creating a function that computes the roots of a quadratic equation
function compute_roots(a, b, c)
    if a ~= 0  %making sure that denominator is nonzero and checking the degree of the equation in hand

       % calculating the roots using the general rule of quadratic equation
       root1 = ((-1 * b) + sqrt((b * b) - (4*a*c)))/ (2 * a); 
       root2 = ((-1 * b) - sqrt((b * b) - (4*a*c)))/ (2 * a);
       fprintf("The first root is %f + %f i \n",real(root1), imag(root1));
       fprintf("The second root is %f + %f i \n",real(root2), imag(root2));
          
    else
      z % the equation is linear in this case (a = 0)
       root1 = -1*c/b; 
       root2 = root1;
       fprintf("The equation is linear and there is only one solution: %f \n",root2);
    end
        
   
end
compute_roots(a1, b1, c1);
compute_roots(a2, b2, c2);
compute_roots(a3, b3, c3);


%% Question 2
%{
You are asked to compute e^(-2) using two different series expansions and compare the results with the true value, which is approximately 1.35335×10−11.35335 \times 10^{-1}1.35335×10−1.
Series Approximations:
First Series (for e^(-x)):
e^(-x) = 1 - x + (x^2)/2! - (x^3)/3! + ...
Second Series (for e^x and then e^(-x)):
e^x=1 + x + (x^2)/2! + (x^3)/3! ...
Then evaluate e^(-x) = 1/(e^x)​.
Steps:
Use 15 terms from each series to compute e−2e^{-2}e−2.
For each series, calculate the true relative error (EtE_tEt​) and the approximate relative error (EaE_aEa​) after each term is added.
Present the results in a table, where each row corresponds to the number of terms used. In the table, include the following columns:
Number of terms (1 through 15)
Approximate value of e−2e^{-2}e−2 for each series
True relative error EtE_tEt​ (in percent)
Approximate relative error EaE_aEa​ (in percent)

}%
%% using the first approach
z = 2;
true_value = 1.35335 * power(10, -1);
num_terms = 1;
% initializing the values of the errors
approx_error = 0;
true_error = 0;
expon_z = 1;
while num_terms <= 15 %iterate for 15 terms
    disp(num_terms)
    % adding new extra term
    new_expon_z = expon_z + ((power(-1, num_terms)*power(z, num_terms))/factorial(num_terms)) 

    true_error = ((true_value - new_expon_z)/true_value)*100 
    approx_error = ((new_expon_z - expon_z)/new_expon_z)*100

    num_terms = num_terms + 1;
    expon_z = new_expon_z;
end
%% Question 2: using the second approach
z = 2;
true_value = 1.35335 * power(10, -1);
num_terms = 1;
approx_error = 0;
true_error = 0;
expon_z = 1; %initializing the value of e^z
expon_neg_z = 1 %initializing the value of e^-z
while num_terms <= 15
    disp(num_terms)
    % adding new extra term
    new_expon_z = expon_z + (power(z, num_terms)/factorial(num_terms)) %adding a new term to e^z   
    expon_neg_z = 1/new_expon_z

    true_error = ((true_value - expon_neg_z)/true_value)*100
    approx_error = ((expon_neg_z - (1/expon_z))/expon_neg_z)*100

    num_terms = num_terms + 1;
    expon_z = new_expon_z;
end
%% Question 5
%{
You are asked to determine how many terms of the Maclaurin series are required to approximate cos⁡(x)\cos(x)cos(x) to 6 significant figures. The Maclaurin series for cos⁡(x)\cos(x)cos(x) is given by:
cos⁡(x) = 1 - (x^2)/2! + (x^4)/4! - (x^6)/6! + (x^8)/8! ...
Steps:
Calculate the approximation for cos⁡(x)\cos(x)cos(x) using x=0.4πx = 0.4\pix=0.4π.
Write a program to determine how many terms are needed to achieve the desired accuracy of 6 significant figures.
Hint:
Recall the formula for the stopping condition based on the relative approximate error ϵs\epsilon_sϵs​, which is often used when aiming for a specific number of significant digits. This concept was discussed in class and is also covered in the textbook.
}%
x = 0.4 * pi;
num_terms = 1;
num_significant_figs = 6;
error_criterion = 0.5 * power(10, 2-num_significant_figs);
approx_error = 10000;
cos_x = 1;
while approx_error > error_criterion
    new_cos_x = cos_x + ((power(-1, num_terms)*power(x, 2*num_terms))/factorial(2*num_terms))
    %calculates the absolute value of
    % the error to make sure that its magnitude will be less than the criterion error value
    approx_error = abs(((new_cos_x - cos_x)/new_cos_x))*100
    num_terms = num_terms + 1
    cos_x = new_cos_x;
end
disp(num_terms)
