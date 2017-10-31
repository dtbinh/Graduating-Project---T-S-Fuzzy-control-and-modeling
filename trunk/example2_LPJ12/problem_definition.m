function [xi, A, x1, x2, h, n] = problem_definition()
% Let us consider a T-S fuzzy model with vertices
% A1 = [-2 4; -1 -2]; A2 = [-2 4; -(1+lambda) -2]. With lambda = 20;
% h1(z(t)) = alpha1(z(t)) = (1 + sin(x1(t)))/2;
% h2(z(t)) = alpha2(z(t)) = 1 - alpha1(z(t));
% C1 = {x(t) E R^n ||xi(t)| <= pi/2, i = 1, 2}.

xi = [-pi/2; pi/2];
lambda = 20;
A(:, :, 1) = [-2 4; -1 -2];
A(:, :, 2) = [-2 4; -(1+lambda) -2];
syms x1 x2;
h(1) = (1 + sin(x1))/2;
h(2) = 1 - h(1);

% number of state variables
n = 2;

end

