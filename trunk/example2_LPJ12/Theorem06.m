function LMIs = Theorem06(A, Pi, n, xi, h, x1);

for i = 1:n
    dh(i) = diff(h(i));
end

points = 100;
x_range = linspace(xi(1), xi(2), points);
syms y;
for i = 1:points
    for j = 1:n
        s = vpasolve([y == dh(j), x1 == x_range(i)]);
        phi_i(i,j) = abs(double(s.y));
    end
end

LMIs = [];

% Pi = Pi' > 0
for i = 1:n
    LMIs = [LMIs, Pi{i} > 0];
end

% Pi + X >= 0
P_phi = 0;
X = sdpvar(n, n, 'symmetric');
for i = 1:n
    LMIs = [LMIs, Pi{i} + X >= 0];
    phi(i) = max(phi_i(:,i));
    P_phi = P_phi + phi(i) * (Pi{i} + X);
end

% P_phi + (1/2) * (Ai'*Pj + Pj*Ai +Aj'*Pi + Pi*Aj) <= 0
for i = 1:n
    for j = 1:n
        LMIs = [LMIs, ]
    end
end

end
