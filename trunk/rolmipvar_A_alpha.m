function poly_A = rolmipvar_A_alpha(Ai, n_alpha, n_theta, n_gamma)

% This function returns the rolmipvar construction of A(alpha, theta,
% gamma). It is known that A depends only of alpha, so de degrees of theta
% and gamma are both zero. It is also known that A has n vertices, each
% one of that associated to a different alpha that when are added together
% resuts on A.

theta = zeros(1, n_theta);
gamma = zeros(1, n_gamma);

for i = 1:n_alpha
    alpha = zeros(1, n_alpha);
    alpha(i) = 1;
    A{i} = {alpha, theta, gamma, Ai(:, :, i)};
end

poly_A = rolmipvar(A,'A', [n_alpha n_theta n_gamma], [1 0 0]);
end