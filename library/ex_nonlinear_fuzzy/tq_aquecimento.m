function output = tq_aquecimento
%close all; clc;

% Definição das constantes do modelo
R = 0.3;  % h/m2  
A = 2; % m2
cp = 0.75;  % kJ/(kg . K)
rho = 1000; % kg/m3 
U = 150;  % kJ/(m2 . s . K)
param = [U A rho cp R];

% Tempo de simulação
ts = 0.0:0.01:60.0; % h
%ts = [0 40];

% Entradas (disturbios) nos instantes %t1 = 2; t2 = 30;
Fe0 = 10; % m3/h  
Te0 = 530; % K  
Th0 = 540;  % K
inputs = [Fe0 Te0 Th0];

% Regime permanente
h0 = (Fe0*R)^2; 
%h0 = (5/A); 
T0 = (Fe0*Te0*rho*cp+U*Th0*A)/(Fe0*rho*cp+U*A);
%T0 = Th(1);

sparam = [U A rho cp R Fe0 Te0 Th0 h0 T0];


% Simulação dinamica não-linear
odesim = 1;
options = odeset('RelTol',1e-4,'AbsTol',1e-9);
[t,x] = ode45('dxdt_tq_aquecimento',ts,[h0 T0],options,[param inputs odesim]);
output.t = t;
output.h = x(:,1);
output.T = x(:,2);

% Simulação dinamica fuzzy TS

% parameters: fuzzy linear matrices
global Ag Bg

a=rand(2); a=a-(max(eig(a))+0.5)*eye(2);eig(a);
Ag{1} = a; Ag{2} = a*2;
Ag{3} = a*1/2; Ag{4} = a*5;
Bg{1} = rand(2,3); Bg{2} = rand(2,3);
Bg{3} = rand(2,3); Bg{4} = rand(2,3);

odesim = 1;
options = odeset('RelTol',1e-4,'AbsTol',1e-9);
[T,X] = ode45('dxdt_tq_aquecimento_fuzzy',ts,[h0 T0],options,[inputs odesim]);
output.t = T;
output.h = X(:,1);
output.T = X(:,2);


% Visualização da simulação
figure(1);
plot(t,x(:,1)); title ('Tanque de aquecimento' );
xlabel('Tempo (h)'); ylabel ('Altura (m)'); grid;
figure(2);
plot(t,x(:,2)); title ('Tanque de aquecimento' );
xlabel('Tempo (h)'); ylabel ('Temperatura (K)' );
grid;  %axis([ts 8.5 11])




