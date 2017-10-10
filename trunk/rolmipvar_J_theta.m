function poly_J = rolmipvar_J_theta(Ji, n_alpha, n_theta, n_gamma)

% This function returns the rolmipvar construction of J(alpha, theta,
% gamma). It is known that J only depends of theta, so de degrees of alpha
% and gamma are both zero. It is also known that J has n_theta vertices, 
% each one of that associated to a different theta that when are added
% together resuts on J.

alpha = zeros(1, n_alpha);
gamma = zeros(1, n_gamma);

for i = 1:n_theta
    theta = zeros(1, n_theta);
    theta(i) = 1;
    J{i} = {alpha, theta, gamma, Ji(:, :, i)};
end

poly_J = rolmipvar(J,'J', [n_alpha n_theta n_gamma], [0 1 0]);
end

