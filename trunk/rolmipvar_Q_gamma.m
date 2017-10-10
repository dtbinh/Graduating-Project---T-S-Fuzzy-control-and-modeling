function poly_Q = rolmipvar_Q_gamma(Qi, n_alpha, n_theta, n_gamma)

% This function returns the rolmipvar construction of Q(alpha, theta,
% gamma). It is known that Q depends only of gamma, so de degrees of alpha
% andtheta are both zero. It is also known that Q has n_gamma vertices,
% each one of that associated to a different gamma that when are added
%together resuts on Q.

alpha = zeros(1, n_alpha);
theta = zeros(1, n_theta);

for i = 1:n_gamma
    gamma = zeros(1, n_gamma);
    gamma(i) = 1;
    Q{i} = {alpha, theta, gamma, Qi{i}};
end

poly_Q = rolmipvar(Q,'Q', [n_alpha n_theta n_gamma], [0 0 1]);
end