function LMIs = Theorem06(LMIs, A, Pi, n, phi)

% Pi = Pi' > 0
for i = 1:n
    LMIs = [LMIs, Pi{i} > 0];
end

% Pi + X >= 0
P_phi = 0;
X = sdpvar(n, n, 'symmetric');
for i = 1:n
    LMIs = [LMIs, Pi{i} + X >= 0];
    P_phi = P_phi + phi(i) * (Pi{i} + X);
end

% P_phi + (1/2) * (Ai'*Pj + Pj*Ai +Aj'*Pi + Pi*Aj) <= 0
for i = 1:n
    for j = 1:n
        LMIs = [LMIs, P_phi + (1/2) * (A(:, :, i)'*Pi{j} + Pi{j}*A(:, :, i) +A(:, :, j)'*Pi{i} + Pi{i}*A(:, :, j)) <= 0];
    end
end

end