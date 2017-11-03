function A = linear_model(A_fuzzy, h, x1)

% nonlinear model
A_ = h(1)*A_fuzzy(:, :, 1) + h(2)*A_fuzzy(:,:,2);

%linearization
A_lin = taylor(A_, x1);
syms y;
s = vpasolve([y == A_lin(2, 1), x1 == 0]);
A = [A_lin(1, 1) A_lin(1, 2); double(s.y) A_lin(2,2)];
A = eval(A);

end

