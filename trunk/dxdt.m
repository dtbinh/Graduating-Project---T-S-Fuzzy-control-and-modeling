function DXDT = dxdt(t,x) 

% states
% x1 = Pf    \in [0; 25000]
% x2 = Qf    \in [-70000; 5000]
% x3 = delta \in [-0.02; 0.1]

[wf, V, Eref, m, n, Ro, Pref, Qref] = parameters;

dx1dt = -wf*x(1) + (wf*(V*(Eref-n*(x(1)-Pref))*cos(x(3))-V^2))/Ro; 
dx2dt = -wf*x(2) - (wf*V*(Eref-n*(x(1)-Pref))*sin(x(3)))/Ro; 
dx3dt = m*(x(2)-Qref); 

DXDT = [dx1dt dx2dt dx3dt]';


