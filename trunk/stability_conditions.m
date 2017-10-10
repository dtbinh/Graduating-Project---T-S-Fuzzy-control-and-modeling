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

A = vertices(z_lim);
% n_alpha is the number of vertices of A(alpha)
[~, ~, n_alpha] = size(A);

J = jacobians_vertices();
%n_theta is the number of vertices of J(theta)
[~, ~, n_theta] = size(J);
n_theta = 2;

x_k = polyhedral_set();
[~, n_gamma] = size(x_k);

poly_A = rolmipvar_A_alpha(A, n_alpha, n_theta, n_gamma);

poly_J = rolmipvar_J_theta(J, n_alpha, n_theta, n_gamma);

for i = 1:n_alpha
    Pi{i} = sdpvar(3,3,'symmetric');
end
poly_P = rolmipvar_P_alpha(Pi, n_alpha, n_theta, n_gamma);

for k = 1:n_gamma
    for i = 1:n_alpha
        if i == 1
            Qi{k} = Pi{i}*x_k(:, n_gamma);
        else
            Qi{k} = horzcat(Qi{k}, Pi{i}*x_k(:, n_gamma));
        end
    end
end

poly_Q = rolmipvar_Q_gamma(Qi, n_alpha, n_theta, n_gamma);

for i = 1:n_alpha
    Xi{i} = sdpvar(3,3,'full');
end
poly_X = rolmipvar_X_alpha(Xi, n_alpha, n_theta, n_gamma);

end

