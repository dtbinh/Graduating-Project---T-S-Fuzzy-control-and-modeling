function Ai = vertices(z_lim)

global X

% nao linearidades
%
% z1 = cos(x3)
% z2 = phi(x3)/x3
% z3 = sen(x3)
% z4 = sen(x3)/x3

[wf, V, Eref, m, n, Ro, Pref, Qref] = parameters;

Ai = zeros(3, 3, 2^4);

i = 0;
for j = 1:2
    for k = 1:2
        for l = 1:2
            for p = 1:2
            i = i + 1;
            %vertices
            Ai(:, :, i) = [-wf-(wf*V*n)/Ro*z_lim{1,j}  0   z_lim{2,k};
                           (wf*V*n)/Ro*z_lim{3,l} -wf -(wf*V*(Eref-n*(X(1)-Pref)))/Ro*z_lim{4,p};
                           0  m   0];
            end
        end
    end
end


