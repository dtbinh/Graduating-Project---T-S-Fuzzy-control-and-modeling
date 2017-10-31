function x_k = polyhedral_set()

% This function consists on representing the k vertices
% that describes the domain of validity of the TS model.
% This representation corresponds to a polyhedral set.

% Physical boundaries of state variables
% x1 = Pf    \in [0; 25000]
% x2 = Qf    \in [-70000; 5000]
% x3 = delta \in [-0.02; 0.1]

% Breakeven point
% X = xe = [1.6252e+04 0 0]'
global X

% As the breakeven point is different of zero, so the it is
% needed to use deviation variables as follows.
% xd = x - xe.

xd1 = [0 - X(1); 25000 - X(1)];
xd2 = [-70000 - X(2); 5000 - X(2)];
xd3 = [-0.02 - X(3); 0.1 - X(3)];

% The number of elements of x_k will corresponds to the
% number of all possible combinations of the boundaries
% given by the minimum and maximum values that is possible
%to each state assumes.

l = 1;
for i = 1:2
    for j = 1:2
        for k = 1:2
            x_k(1, l) = xd1(i);
            x_k(2, l) = xd2(j);
            x_k(3, l) = xd3(k);
            l = l + 1;
        end
    end
end

end

