%clear all; clear; clc;

% Let us consider a T-S fuzzy model with vertices
% A1 = [-2 4; -1 -2]; A2 = [-2 4; -(1+lambda) -2]. With lambda = 20;
% h1(z(t)) = alpha1(z(t)) = (1 + sin(x1(t)))/2;
% h2(z(t)) = alpha2(z(t)) = 1 - alpha1(z(t));
% C1 = {x(t) E R^n ||xi(t)| <= pi/2, i = 1, 2}.

xi = [-pi/2; pi/2];
lambda = 20;
A_(:, :, 1) = [-2 4; -1 -2];
A_(:, :, 2) = [-2 4; -(1+lambda) -2];
syms x1 x2
h1 = (1 + sin(x1))/2;
h2 = 1 - h1;

x_k = StateVariablesVertices(xi);

% method 4) Fuzzy dynamics TS + Lyapunov method with P(alpha): LMIs (8)-(9)

% LMI (8)
[LMIs_4_1, P, n_alpha, n_theta, n_gamma] = ...
                                fuzzy_TS_dinamic_Lyapunov_with_P_alpha(...
                                                A_, h1, h2, x1, x2, x_k);
% LMI (9)                                            
LMIs_4_2 = LargestInvariantSetContainedInPolytope(LMIs_4_1, x_k, P);

% sol = solvesdp(LMIs_4_2,[],sdpsettings('verbose',0,'solver','sedumi'));
% %retrieving the minimal primal residual
% p = min(checkset(LMIs_4_2));
% display(p);
% 
% %capturing the solutions (if ones exist)
% maxViolation = 1e-7; %minimization problem
% %maxViolation = 0;   %factibility problem
% if p  > -maxViolation %adopted precision for the minimum primal residual
%     msgbox 'Stable (method 4)';
%     output.P = double(P);
%     P_n = verticesP(output.P, n_alpha, n_theta, n_gamma);
%     level_curve(P_n, 'r');
% else
%     msgbox 'Not stable (method 4) <enter>';
% end

% method 4) + set enlargement
[LMIs_4_3, crit] = EnlargementOfLargestInvariantSet(LMIs_4_2, P);

% solvesdp(LMIs_4_3,crit,sdpsettings('solver','sedumi','verbose',0));
% [p,d]=checkset(LMIs_4_3);
% pmin = min(checkset(LMIs_4_3));
% display(pmin)
% if sum(p > -maxViolation)
% 	msgbox 'Stable  (method 4 + set enlargement')'
%     output.P = double(P);
%     P_n = verticesP(output.P, n_alpha, n_theta, n_gamma);
%     level_curve(P_n, 'g');
% else
%     msgbox 'Not stable (method 4 + set enlargement)'
% end

% theorem 2
[LMIs, T] = theorem2(xi, P, LMIs_4_3);
omega1 = 1;
omega2 = 1;
crit = omega1*trace(T)+omega2*gamma;
%crit = [];

solvesdp(LMIs,crit,sdpsettings('solver','sedumi','verbose',0));
[p,d]=checkset(LMIs);
pmin = min(checkset(LMIs));
display(pmin)

tol = 1e-7;
if sum(p < -tol)
    msgbox 'instavel'
end
P=double(P);
W = inv(P);
gamma=double(gamma);

projellisa(P*gamma,'b','-')

volume = sqrt(det(W/gamma))