clear all; clear; clc;

lambda = 20;
[xi, A, x1, x2, h, n] = problem_definition(lambda);

x_k = StateVariablesVertices(xi);

% method 2) Fuzzy dynamics TS + Lyapunov method with P constant: LMIs (9)
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
if pmin > -maxViolation
	msgbox 'Stable  (method 2')'
    output.P = double(poly_P);
    P_n = {};
    P_n{1} = output.P;
    level_curve(P_n, 1, 'r');
else
    msgbox 'Not stable (method 2)'
end
