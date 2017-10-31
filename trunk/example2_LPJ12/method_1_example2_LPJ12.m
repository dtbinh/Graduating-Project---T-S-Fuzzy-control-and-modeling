% Let us consider a T-S fuzzy model with vertices
% A1 = [-2 4; -1 -2]; A2 = [-2 4; -(1+lambda) -2]. With lambda = 20;
% h1(z(t)) = alpha1(z(t)) = (1 + sin(x1(t)))/2;
% h2(z(t)) = alpha2(z(t)) = 1 - alpha1(z(t));
% C1 = {x(t) E R^n ||xi(t)| <= pi/2, i = 1, 2}.

% number of state variables
n = 2;
xi = [-pi/2; pi/2];
x_k = StateVariablesVertices(xi);

lambda = 20;
A_fuzzy(:, :, 1) = [-2 4; -1 -2];
A_fuzzy(:, :, 2) = [-2 4; -(1+lambda) -2];
syms x1 x2;
h(1) = (1 + sin(x1))/2;
h(2) = 1 - h(1);

A_ = h(1)*A_fuzzy(:, :, 1) + h(2)*A_fuzzy(:,:,2);
A_lin = taylor(A_, x1);
syms y;
s = vpasolve([y == A_lin(2, 1), x1 == 0]);
A = [A_lin(1, 1) A_lin(1, 2); double(s.y) A_lin(2,2)];
A = eval(A);

P = sdpvar(n,n,'symmetric');

LMIs = [];
LMIs = [LMIs, P >= 0];
LMIs = [LMIs, A'*P + P*A <= 0];

LMIs = LargestInvariantSetContainedInPolytope(LMIs, x_k, P);

sol = solvesdp(LMIs,[],sdpsettings('verbose',0,'solver','sedumi'));
%retrieving the minimal primal residual
p = min(checkset(LMIs));
display(p);

%capturing the solutions (if ones exist)
maxViolation = 1e-7; %minimization problem
%maxViolation = 0;   %factibility problem
if p  > -maxViolation %adopted precision for the minimum primal residual
    msgbox 'Stable (method 1)';
    output.P = double(P);
    P_n = verticesP(output.P, n_alpha, n_theta, n_gamma);
    level_curve(P_n, 'y');
else
    msgbox 'Not stable (method 1)';
end