function [LMIs, T, gamma] = theorem2(xi, P, LMIs, n_alpha, n_theta, n_gamma, n)

mu = [max(xi) max(xi)];
In = eye(n);
gamma = sdpvar(1,1,'symmetric');

% T
T = rolmipvar(n,n,'T','full',[n_alpha n_theta n_gamma], [0 0 0]);
LMIs = [LMIs, [T P;P P]>=0];

for i = 1:n
        MatQ = [P*mu(i) In(i,:)'; In(i,:) gamma*mu(i)];
        LMIs = [LMIs, MatQ >= 0];
end


end