function poly_P = rolmipvar_P_alpha(Pi, n_alpha, n_theta, n_gamma)

% This function returns the rolmipvar construction of P(alpha, theta,
% gamma). It is known that P depends only of alpha, so de degrees of theta
% and gamma are both zero. It is also known that P has n_gamma vertices,
% each one of that associated to a different alpha that when are added
% together resuts on P.

theta = zeros(1, n_theta);
gamma = zeros(1, n_gamma);

for i = 1:n_alpha
    alpha = zeros(1, n_alpha);
    alpha(i) = 1;
    P{i} = {alpha, theta, gamma, Pi{i}};
end

poly_P = rolmipvar(P,'P', [n_alpha n_theta n_gamma], [1 0 0]);
end
