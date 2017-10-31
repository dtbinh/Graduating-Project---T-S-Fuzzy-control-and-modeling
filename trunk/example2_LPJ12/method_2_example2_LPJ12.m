%clear all; clear; clc;

% Let us consider a T-S fuzzy model with vertices
% A1 = [-2 4; -1 -2]; A2 = [-2 4; -(1+lambda) -2]. With lambda = 20;
% h1(z(t)) = alpha1(z(t)) = (1 + sin(x1(t)))/2;
% h2(z(t)) = alpha2(z(t)) = 1 - alpha1(z(t));
% C1 = {x(t) E R^n ||xi(t)| <= pi/2, i = 1, 2}.

xi = [-pi/2; pi/2];
lambda = 20;
A1 = [-2 4; -1 -2];
A2 = [-2 4; -(1+lambda) -2];
syms x1 x2
h1 = (1 + sin(x1))/2;
h2 = 1 - h1;

% number of state variables
n = 2;
x_k = StateVariablesVertices(xi);

% method 2) Fuzzy dynamics TS + Lyapunov method with P constant: LMIs (9)
A = [A1 A2];
poly_A = rolmipvar(A,'A',2,1);
P = sdpvar(n,n,'symmetric');
poly_P = rolmipvar(P,'P',2,0);

LMIs = [];
LMIs = [LMIs, poly_A'*poly_P + poly_P*poly_A <= 0];
LMIs = LargestInvariantSetContainedInPolytope(LMIs, x_k, poly_P);
[LMIs, crit] = EnlargementOfLargestInvariantSet(LMIs, poly_P);

solvesdp(LMIs, crit, sdpsettings('solver', 'sedumi', 'verbose', 0));
[p,d]=checkset(LMIs);
pmin = min(checkset(LMIs));
display(pmin)
maxViolation = 1e-7; %minimization problem
if sum(p > -maxViolation)
	msgbox 'Stable  (method 4 + set enlargement')'
    output.P = double(poly_P);
    P_n{1} = output.P;
    level_curve(P_n, 'y');
else
    msgbox 'Not stable (method 4 + set enlargement)'
end
