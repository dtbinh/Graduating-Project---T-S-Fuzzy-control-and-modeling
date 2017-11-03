function [dh, phi_max, phi_min] = determination_of_phi_range_and_diff_h(n, h, A, xi, x1, x2)

A_nonlinear = 0;
for i = 1:n
    A_nonlinear = A_nonlinear + h(i)*A(:, :, i);
end

for i = 1:n
    dh(i) = diff(h(i)) * (A_nonlinear(1, :) * [x1; x2]);
end

points = 10;
x_range = linspace(xi(1), xi(2), points);
syms y;
for i = 1:points
    for j = 1:points
        for k = 1:n
            s = vpasolve([y == dh(k), x1 == x_range(i), x2 == x_range(j)]);
            phi_i_j(i, j,k) = abs(double(s.y));
        end
        for k = 1:n
            phi_i_max(i, k) = max(phi_i_j(i, :,k));
            phi_i_min(i, k) = min(phi_i_j(i, :,k));
        end
    end
end

for i = 1:n
    phi_min(i) = min(phi_i_min(:, i));
    phi_max(i) = max(phi_i_max(:, i));
%     phi_max(i) = 9.4248;
%     phi_min(i) = 0;
end



end

