function [DXDT_FUZZY ] = dxdt_fuzzy(t, x)


% states
% x1 = Pf    \in [0; 25000]
% x2 = Qf    \in [-70000; 5000]
% x3 = delta \in [-0.02; 0.1]

% nao linearidades
%
% z1 = cos(x(3))
% z2 = (-wf*(1 + V*n*cox(x3)/Ro))
%       + (wf*V*(Eref + n*(-Pf0 + Pref))*(cos(x3)/x3)/Ro)*X(1)
%       - wf*(X(1)+V^2/Ro)
% z3 = sen(x3)
% z4 = sen(x3)/x3

% limites das nao-linearidades
global z_lim

% numero de nao-linearidades
N = 4;

% graus de petinencia
Mij = membership_degrees(x, N);

% funcoes de pertinencia
alpha = membership_functions(Mij, N);

%vertices
Ai = vertices(z_lim);

% modelo fuzzy
dxdt = zeros(3, 1);

for i = 1:2^N
    dxdt = dxdt + alpha(i)*(Ai(:, :, i)*x);
end

DXDT_FUZZY = [dxdt(1,:) dxdt(2,:) dxdt(3,:)]';


