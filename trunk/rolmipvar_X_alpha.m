function poly_X = rolmipvar_X_alpha(Xi, n_alpha, n_theta, n_gamma)

% This function returns the rolmipvar construction of X(alpha, theta,
% gamma). It is known that X depends only of alpha, so de degrees of theta
% and gamma are both zero. It is also known that X has n_alpha vertices,
% each one of that associated to a different alpha that when are added
% together resuts on X.

theta = zeros(1, n_theta);
gamma = zeros(1, n_gamma);

for i = 1:n_alpha
    alpha = zeros(1, n_alpha);
    alpha(i) = 1;
    X{i} = {alpha, theta, gamma, Xi{i}};
end

poly_X = rolmipvar(X, 'X', [n_alpha n_theta n_gamma], [1 0 0]);
end
