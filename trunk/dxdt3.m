%function DXDT = dxdt2() 
function DXDT = dxdt3(t,x) 

global X

% states
% x1 = Pf    \in [0; 25000]
% x2 = Qf    \in [-70000; 5000]
% x3 = delta \in [-0.02; 0.1]

%syms x3 x2 x1,
x1 = x(1);
x2 = x(2);
x3 = x(3);

[wf, V, Eref, m, n, Ro, Pref, Qref] = parameters;

x1d = x1 + X(1);
x2d = x2 + X(2);
x3d = x3 + X(3);

dx1dt = -wf*x1d + (wf*(V*(Eref-n*(x1d-Pref))*cos(x3d)-V^2))/Ro;
dx2dt = -wf*x2d - (wf*V*(Eref-n*(x1d-Pref))*sin(x3d))/Ro;
dx3dt = m*(x2d-Qref);

DXDT = [dx1dt dx2dt dx3dt]';