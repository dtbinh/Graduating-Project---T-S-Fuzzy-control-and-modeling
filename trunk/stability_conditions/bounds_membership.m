function z_lim = bounds_membership(X)

[wf, V, Eref, m, n, Ro, Pref, Qref] = parameters;

% variacoes das nao linearidades
% em funcao de (delta = x3)
N = 100;
x3_range = linspace(-0.02, 0.1, N);
x1_range = linspace(-X(1), 25000-X(1), N); % desprezou-se o ponto x1 = 0 para evitar indeterminacao

% faixas de atuacao das nao linearidades
for i = 1:N
    z1(i) = cos(x3_range(i));
    z2(i) = ((wf*(V*(Eref-n*(X(1)-Pref))*cos(x3_range(i))-V^2))/Ro -wf*X(1))/x3_range(i);
    z3(i) = sin(x3_range(i));
    z4(i) = sin(x3_range(i))/x3_range(i);
end

% z1 = cos(x3)
z_lim{1,1} = max(z1);
z_lim{1,2} = min(z1);

%z2 = phi(x3)/x3
z_lim{2,1} = max(z2);
z_lim{2,2} = min(z2);

%z3 = sin(x3)
z_lim{3,1} = max(z3);
z_lim{3,2} = min(z3);

%z4 = sin(x3)/x3
z_lim{4,1} = max(z4);
z_lim{4,2} = min(z4);