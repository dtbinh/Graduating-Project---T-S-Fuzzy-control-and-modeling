function x_k = polyhedral_set()

% This function consists on representing the k vertices
% that describes the domain of validity of the TS model.
% This representation corresponds to a polyhedral set.

% Breakeven point
% x1 = Pf    \in [0; 25000]
% x2 = Qf    \in [-70000; 5000]
% x3 = delta \in [-0.02; 0.1]

x1 = [0; 25000];
x2 = [-70000; 5000];
x3 = [-0.02; 0.1];

% The number of elements of x_k will corresponds to the
% number of all possible combinations of the boundaries
% given by the minimum and maximum values that is possible
%to each state assumes.

l = 1;
for i = 1:2
    for j = 1:2
        for k = 1:2
            x_k(1, l) = x1(i);
            x_k(2, l) = x2(j);
            x_k(3, l) = x3(k);
            l = l + 1;
        end
    end
end

end

