function X = dxdt_fsolve(X0)

if nargin == 0
    N = 0; %0,3,5
    X0 = [1.6252e+004 0 N*pi];
end

X = fsolve(@dxdt,X0);

function DXDT = dxdt(x) 

% states
% x1 = Pf    \in [0; 25000]
% x2 = Qf    \in [-70000; 5000]
% x3 = delta \in [-0.02; 0.1]

wf = 2*pi*60;
V=311;
Eref=311;
m=3.77/22000;
n=20/22000;
Ro=0.1;
Pref=22000;
Qref=0;

dx1dt = -wf*x(1) + (wf*(V*(Eref-n*(x(1)-Pref))*cos(x(3))-V^2))/Ro; 
dx2dt = -wf*x(2) - (wf*V*(Eref-n*(x(1)-Pref))*sin(x(3)))/Ro; 
dx3dt = m*(x(2)-Qref); 

DXDT = [dx1dt dx2dt dx3dt]';