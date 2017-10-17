function stability_conditions()

global z_lim X;

% at this function the following LMIs will be solved, if there is
% an existing solution, to verify if the nonlinear system is stable.
% If there is not a solution (an exiting matrix P) so the system is
% not stable.
%
% A(alpha)' P(alpha) + P(alpha) A(alpha) + Q(gamma) J(theta) A(alpha) + ...
% ... X(alpha) 1' J(theta) A(alpha) + A(alpha)' J(theta)' 1 X(alpha)' < 0
%
% Where Q = sum(gamma^k *[P_1*x^k ... P_N*x^k])
% the vertices of X(alpha) and P(alpha) are the variables

% Breakeven point - initial conditions
X = dxdt_fsolve;
z_lim = bounds_membership(X);

A_alpha = vertices(z_lim);
% n_alpha is the number of vertices of A(alpha)
[~, ~, n_alpha] = size(A_alpha);

J_theta = jacobians_vertices();
%n_theta is the number of vertices of J(theta)
[~, ~, n_theta] = size(J_theta);

x_k = polyhedral_set();
[~, n_gamma] = size(x_k);

A = rolmipvar_Matrix(A_alpha, 'A', n_alpha, n_theta, n_gamma, [1 0 0]);

J = rolmipvar_Matrix(J_theta, 'J', n_alpha, n_theta, n_gamma, [0 1 0]);

theta = zeros(1, n_theta);
gamma = zeros(1, n_gamma);
for i = 1:n_alpha
    alpha = zeros(1, n_alpha);
    alpha(i) = 1;
    Pi_{i} = sdpvar(3,3,'symmetric');
    Pi{i} = {alpha, theta, gamma, Pi_{i}};
end
poly_P = rolmipvar(Pi,'P', [n_alpha n_theta n_gamma], [1 0 0]);

alpha = zeros(1, n_alpha);
theta = zeros(1, n_theta);
for k = 1:n_gamma
    for i = 1:n_alpha
        if i == 1
            Qi_{k} = Pi_{i}*x_k(:, n_gamma);
        else
            Qi_{k} = horzcat(Qi_{k}, Pi_{i}*x_k(:, n_gamma));
        end
    end
    gamma = zeros(1, n_gamma);
    gamma(k) = 1;
    Qi{k} = {alpha, theta, gamma, Qi_{k}};
end
Q = rolmipvar(Qi,'Q', [n_alpha n_theta n_gamma], [0 0 1]);

theta = zeros(1, n_theta);
gamma = zeros(1, n_gamma);
for i = 1:n_alpha
    alpha = zeros(1, n_alpha);
    alpha(i) = 1;
    Xi{i} = {alpha, theta, gamma, sdpvar(3,3,'full')};
end
X_ = rolmipvar(Xi, 'X', [n_alpha n_theta n_gamma], [1 0 0]);

LMIs = [ [A'*poly_P + poly_P*A + Q*J*A + X_* ones(16,3)'*J*A + A'*J'*ones(16,3)*X_']< 0];

poly = polytope(x_k');
[H,K] = double(poly);
[q, ~] = size(H);
b_k = [];
for k = 1:q
    b_k = [b_k (H(k, :)/K(k))'];
end
[r, c] = size(b_k');
LMIs = [LMIs, [ones(max(r, c), max(r, c)) b_k' ;b_k poly_P] >=0];

sol = solvesdp(LMIs,[],sdpsettings('verbose',0,'solver','sedumi'));

% retrieving the minimal primal residual
p=min(checkset(LMIs));
display(p);

%capturing the solutions (if ones exist)
maxViolation = 1e-7; %minimization problem
%maxViolation = 0;   %factibility problem
if p  > -maxViolation %adopted precision for the minimum primal residual
    output.P = double(poly_P);
end

end

