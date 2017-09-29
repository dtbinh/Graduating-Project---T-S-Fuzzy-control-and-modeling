function dx = dxdt_tq_aquecimento(t,x,flag,par)

% parameters
U = par(1);
A = par(2);
Ro = par(3);
Cp  = par(4);
R = par(5);

%inputs
Fe  = par(6);
Te = par(7);
Th  = par(8);

odesim = par(9);

if odesim == 1
    Fei = f_input(t,Fe);
    Tei = f_input(t,Te);
    Thi = f_input(t,Th);
else
    Fei  = par(6);
    Tei = par(7);
    Thi  = par(8);
end

%function
dh = (Fei-(sqrt(x(1))/R))/A;
dT = (1/x(1))*((Fei*Tei/A)+(U*Thi/(Ro*Cp)) - ...
((Fei/A)+(U/(Ro*Cp)))*x(2));

dx  = [dh; dT];


%___________________________________
function fe = f_input(t,x)
% disturbio de entrada
% x = [x1 x2 t1 t2]
if t < 2 || t > 30
    fe = x;
else
    fe = x*1.2;
end

