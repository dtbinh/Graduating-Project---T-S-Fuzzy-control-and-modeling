%clear all; clear; clc;

% Let us consider a T-S fuzzy model with vertices
% A1 = [-2 4; -1 -2]; A2 = [-2 4; -(1+lambda) -2]. With lambda = 20;
% h1(z(t)) = alpha1(z(t)) = (1 + sin(x1(t)))/2;
% h2(z(t)) = alpha2(z(t)) = 1 - alpha1(z(t));
% C1 = {x(t) E R^n ||xi(t)| <= pi/2, i = 1, 2}.

xi = [-pi/2; pi/2];
lambda = 20;
Ai(:, :, 1) = [-2 4; -1 -2];
Ai(:, :, 2) = [-2 4; -(1+lambda) -2];
syms x1 x2
h1 = (1 + sin(x1))/2;
h2 = 1 - h1;

% number of state variables
n = 2;
x_k = StateVariablesVertices(xi);

% method 3) Fuzzy dynamics TS + Lyapunov method with P(alpha)
%           + Theorem06, MPS09 + LMIs (9)
A = [Ai(:, :, 1) Ai(:, :, 2)];
poly_A = rolmipvar(A,'A',2,1);
for i = 1:n
    alpha = zeros(1, n);
    alpha(i) = 1;
    Pi{i} = sdpvar(n, n, 'symmetric');
    P_{i} = {alpha, Pi{i}};
end
poly_P = rolmipvar(P_,'P', n, 1);

LMIs = [];
LMIs = [LMIs, poly_A'*poly_P + poly_P*poly_A <= 0];
LMIs = LargestInvariantSetContainedInPolytope(LMIs, x_k, poly_P);

sol = solvesdp(LMIs,[],sdpsettings('verbose',0,'solver','sedumi'));
%retrieving the minimal primal residual
p = min(checkset(LMIs));
display(p);

%capturing the solutions (if ones exist)
maxViolation = 1e-7; %minimization problem
%maxViolation = 0;   %factibility problem
if p  > -maxViolation %adopted precision for the minimum primal residual
    msgbox 'Stable (method 4)';
    output.P = double(poly_P);
    for i = 1:n
        alpha = zeros(1, n);
        alpha(i) = 1;
        P_vertices{i} = output.P(alpha);
    end
    level_curve(P_vertices, 'b');
else
    msgbox 'Not stable (method 4)';
end