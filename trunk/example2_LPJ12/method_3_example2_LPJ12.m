clear all; clear; clc;

[xi, A, x1, x2, h, n] = problem_definition();

x_k = StateVariablesVertices(xi);

% method 3) Fuzzy dynamics TS + Lyapunov method with P(alpha)
%           + Theorem06, MPS09 + LMIs (9)
poly_A = rolmipvar(A,'A',2,1);
Pi = {};
P_ = {};
for i = 1:n
    alpha = zeros(1, n);
    alpha(i) = 1;
    Pi{i} = sdpvar(n, n, 'symmetric');
    P_{i} = {alpha, Pi{i}};
end
poly_P = rolmipvar(P_,'P', n, 1);

LMIs = [];
LMIs = [LMIs, poly_A'*poly_P + poly_P*poly_A <= 0];
LMIs = [LMIs, Theorem06(A, Pi, n, xi, h, x1)];
LMIs = LargestInvariantSetContainedInPolytope(LMIs, x_k, poly_P);
[LMIs, crit] = EnlargementOfLargestInvariantSet(LMIs, poly_P);

solvesdp(LMIs, crit, sdpsettings('solver', 'sedumi', 'verbose', 0));
[p,d]=checkset(LMIs);
pmin = min(checkset(LMIs));
display(pmin)
maxViolation = 1e-7; %minimization problem
if sum(p > -maxViolation)
	msgbox 'Stable  (method 3)'
    output.P = double(poly_P);
    P_n = {};
    for i = 1:n
        alpha = zeros(1, n);
        alpha(i) = 1;
        P_n{i} = output.P(alpha);
    end
    level_curve(P_n, 1, 'm');
else
    msgbox 'Not stable (method 3)'
end