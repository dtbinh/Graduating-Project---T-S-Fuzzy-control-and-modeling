function LMIs = Theorem06(Ai, Pi, n)

X = rolmipvar(n,n,'X','full',n, 0);
Xi = X(0);

LMIs = [];

% Pi = Pi' > 0
for i = 1:n
    LMIs = [LMIs, Pi{i} > 0];
end

% Pi = X > 0
for i = 1:n
    LMIs = [LMIs, Pi{i} + Xi > 0];
end

for i = 1:n
    P_phi_i(i) = Pi{i} + Xi;
    alpha = zeros(1, n);
    gamma = zeros(1, n);
    gamma(i) = 1;
    P_{i} = {alpha, gamma, P_phi_i{i}};
end
P_phi = rolmipvar(P_, 'P_phi', [n n], [0 1]);

for i = 1: n
    alpha_i = zeros(1, n);
    alpha_i(i) = 1;
    for j = 1:n
        alpha_j = zeros(1, n);
        alpha_j(j) = 1;
        
    end
end

end

