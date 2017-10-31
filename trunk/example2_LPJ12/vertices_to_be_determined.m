function [ P, Q, X] = vertices_to_be_determined(n_alpha, n_theta, n_gamma,...
                                                x_k, states_number)

theta = zeros(1, n_theta);
gamma = zeros(1, n_gamma);
for i = 1:n_alpha
    alpha = zeros(1, n_alpha);
    alpha(i) = 1;
    Pi{i} = sdpvar(states_number, states_number, 'symmetric');
    P_{i} = {alpha, theta, gamma, Pi{i}};
end
P = rolmipvar(P_,'P', [n_alpha n_theta n_gamma], [1 0 0]);

alpha = zeros(1, n_alpha);
theta = zeros(1, n_theta);
for k = 1:n_gamma
    for i = 1:n_alpha
        if i == 1
            Qi{k} = Pi{i}*x_k(:, n_gamma);
        else
            Qi{k} = horzcat(Qi{k}, Pi{i}*x_k(:, n_gamma));
        end
    end
    gamma = zeros(1, n_gamma);
    gamma(k) = 1;
    Q_{k} = {alpha, theta, gamma, Qi{k}};
end
Q = rolmipvar(Q_,'Q', [n_alpha n_theta n_gamma], [0 0 1]);

theta = zeros(1, n_theta);
gamma = zeros(1, n_gamma);
for i = 1:n_alpha
    alpha = zeros(1, n_alpha);
    alpha(i) = 1;
    Xi{i} = {alpha, theta, gamma, sdpvar(states_number, states_number, 'full')};
end
X = rolmipvar(Xi, 'X', [n_alpha n_theta n_gamma], [1 0 0]);

end

