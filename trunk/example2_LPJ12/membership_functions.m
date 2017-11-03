lambda = 20;
[xi, A_fuzzy, x1, x2, h, n] = problem_definition(lambda);
x_k = StateVariablesVertices(xi);
A = linear_model(A_fuzzy, h, x1);

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
    P_n = {};
    P_n{1} = output.P;
    level_curve(P_n, 1, 'y');
else
    msgbox 'Not stable (method 1)';
end