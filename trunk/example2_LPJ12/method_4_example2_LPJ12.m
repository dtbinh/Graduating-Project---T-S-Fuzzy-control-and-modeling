clear all; clear; clc;

[xi, A, x1, x2, h, n] = problem_definition(20);

x_k = StateVariablesVertices(xi);

% method 4) Fuzzy dynamics TS + Lyapunov method with P(alpha): LMIs (8)-(9)

% LMI (8)
[LMIs, P, n_alpha, n_theta, n_gamma] = ...
                                fuzzy_TS_dinamic_Lyapunov_with_P_alpha(...
                                                A, h(1), h(2), x1, x2, x_k);

applyTheorem2 = false;

if (applyTheorem2 == false)
    % LMI (9)                                            
    LMIs = LargestInvariantSetContainedInPolytope(LMIs, x_k, P);
    [LMIs, crit] = EnlargementOfLargestInvariantSet(LMIs, P);
    solvesdp(LMIs,crit,sdpsettings('solver','sedumi','verbose',0));
    [p,d]=checkset(LMIs);
    pmin = min(checkset(LMIs));
    display(pmin)
    maxViolation = 1e-7;
    if sum(p > -maxViolation)
        msgbox 'Stable  (method 4 + set enlargement')'
        output.P = double(P);
        P_n = verticesP(output.P, n_alpha, n_theta, n_gamma);
        %level_curve(P_n, 1, 'g');
        level_curve_v(P_n, 1, xi, h, 'g');
    else
        msgbox 'Not stable (method 4 + set enlargement)'
    end
else
    % theorem 2
    [LMIs, T, gamma] = theorem2(xi, P, LMIs, n_alpha, n_theta, n_gamma, n);
    omega1 = 1;
    omega2 = 1;
    %gamma = sdpvar(1,1,'symmetric');
    crit = omega1 * trace(T) + omega2 * gamma;
    crit = crit([0]);
    crit = crit{1};

    solvesdp(LMIs,crit,sdpsettings('solver','sedumi','verbose',0));
    [p,d]=checkset(LMIs);
    pmin = min(checkset(LMIs));
    display(pmin)
    maxViolation = 1e-7;
    if pmin > -maxViolation
        msgbox 'Stable  (method 4 (with theorem 2))'
        output.P = double(P);
        P_n = verticesP(output.P, n_alpha, n_theta, n_gamma);
%         level_curve(P_n, 1/double(gamma), 'r');
        level_curve_v(P_n, 1/double(gamma), xi, h, 'r');
    else
        msgbox 'Not stable (method 4 (with theorem 2))'
    end
end