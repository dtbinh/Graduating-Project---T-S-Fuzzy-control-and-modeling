function Mij = membership_degrees(x, N)
% here x = x(3) = delta.

% ponto de equilibrio e limites das nao-linearidades
global X z_lim

[wf, V, Eref, ~, n, Ro, Pref, ~] = parameters();

% graus de petinencia

%Mij=zeros(2,N);
%z  = zeros(1, N);
for j = 1:N
    if j == 1
        z(j) = cos(x(3));
    elseif j == 2
        z(j) = ((wf*(V*(Eref-n*(X(1)-Pref))*cos(x(3))-V^2))/Ro -wf*X(1))/x(3);
    elseif j == 3
        z(j) = sin(x(3));
    elseif j == 4
        z(j) = sin(x(3))/x(3);
    end
    Mij(1,j) = (z_lim{j,1} - z(j))/(z_lim{j,1} - z_lim{j,2});
    Mij(2,j) = (z(j) - z_lim{j,2})/(z_lim{j,1} - z_lim{j,2});
end

end

