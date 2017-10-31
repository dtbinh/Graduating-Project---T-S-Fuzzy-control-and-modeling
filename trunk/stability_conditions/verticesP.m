function P_vertices = verticesP(P, n_alpha, n_theta, n_gamma)

theta = zeros(1, n_theta);
gamma = zeros(1, n_gamma);
for i = 1:n_alpha
    alpha = zeros(1, n_alpha);
    alpha(i) = 1;
    P_vertices{i} = P(alpha, theta, gamma);
end

end

