%function out = main

addpath(genpath('C:\Program Files\MATLAB\YALMIP'))
addpath(genpath('C:\Program Files\MATLAB\SeDuMi_1_3'))
addpath(genpath('C:\Program Files\MATLAB\robust_lmi_parser'))

clear all; close all; clc;

global X z_lim

% Regime permanente
% x1 = Pf    \in [0; 25000]
% x2 = Qf    \in [-70000; 5000]
% x3 = delta \in [-0.02; 0.1]

% constants of the equations
[wf, V, Eref, m, n, Ro, Pref, Qref] = parameters;

% initial conditions
X0 = [10000 1000 0.01];            

Pf0 = X0(1);
Qf0 = X0(2);
delta0 = X0(3);

% Tempo de simulação
ts = [0 1e-1];

%condicoes iniciais ponto de equilibrio
X = dxdt_fsolve;

z_lim = bounds_membership(X);

% solving the nonlinear system
options = odeset('RelTol',1e-4,'AbsTol',1e-9);
[t,x] = ode45('dxdt',ts,[Pf0 Qf0 delta0],options);

out.t = t;
out.naolinear.Pf = x(:,1);
out.naolinear.Qf = x(:,2);
out.naolinear.delta = x(:,3);

% ploting the phase retrait
xend = [];
xend = [xend; x(end,:)];
phase_retrait(x);
disp(xend);

% solving the fuzzy system
[t_fuzzy,x_fuzzy] = ode45('dxdt_fuzzy',ts, [Pf0-X(1) Qf0 delta0], options);
x_fuzzy(:,1) = x_fuzzy(:,1)+X(1);

out.t_fuzzy = t_fuzzy;
out.fuzzy.Pf = x_fuzzy(:,1)+X(1);
out.fuzzy.Qf = x_fuzzy(:,2);
out.fuzzy.delta = x_fuzzy(:,3);

% simulation view
figure;
plot(t,x(:,1),'b',t,x(:,2),'r',t,x(:,3),'g',...
     t_fuzzy,x_fuzzy(:,1),'k',t_fuzzy,x_fuzzy(:,2),'y',t_fuzzy,x_fuzzy(:,3),'c');
xlabel('Tempo'); ylabel ('Pf, Qf, delta'); grid;
legend('Pf', 'Qf', 'delta', 'Pf_fuzzy', 'Qf_fuzzy', 'delta_fuzzy');

figure;
plot(t,x(:,1),t_fuzzy,x_fuzzy(:,1)); title ('P_f' );
figure;
plot(t,x(:,2),t_fuzzy,x_fuzzy(:,2)); title ('Q_f' );
figure;
plot(t,x(:,3),t_fuzzy,x_fuzzy(:,3)); title ('\delta' );