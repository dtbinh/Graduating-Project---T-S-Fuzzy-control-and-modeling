function Ji = jacobians_vertices()

syms x1 x2 x3;

% non-linearities number
N = 4;

% membership degrees
Mij = membership_degrees([x1, x2, x3], N);

% membership functions
alpha = membership_functions(Mij, N);

% jacobians of membership functions
[~, r] = size(alpha);
jacobian = [];
for i = 1: r
    jacobian = [jacobian gradient(alpha(i), [x1, x2, x3])];
end
jacobian = jacobian';

% Permanent regime - physical limits of states
% x1 = Pf    \in [0; 25000]
% x2 = Qf    \in [-70000; 5000]
% x3 = delta \in [-0.02; 0.1]

K = 100;
x3_range = linspace(-0.02, 0.1, K);
[I, J] = size(jacobian);
syms y;
Ji = zeros(I, J, 2);
for i = 1: I
    for j = 1:J
        for k = 1:K
            s = vpasolve([y == jacobian(i, j), x3 == x3_range(k)], [x3, y]);
            J_range(i,j, k) = double(s.y);
        end
        Ji(i, j, 1) = min(J_range(i,j,:));
        Ji(i, j, 2) = max(J_range(i,j,:));
    end
end

end