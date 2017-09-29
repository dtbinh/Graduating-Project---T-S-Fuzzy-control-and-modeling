function dx = dxdt_tq_aquecimento_fuzzy(t,x,flag,par)

% parameters
global Ag Bg


%inputs
Fe  = par(1);
Te = par(2);
Th  = par(3);

odesim = par(4);

if odesim == 1
    Fei = f_input(t,Fe);
    Tei = f_input(t,Te);
    Thi = f_input(t,Th);
else
    Fei  = par(6);
    Tei = par(7);
    Thi  = par(8);
end

u = [Fei Tei Thi]'; %expressar em variaveis de desvio

alphaact = mshipfct(x,20,600); %lmits of x(1) and x(2)
Afuzzy = vert2mu(Ag,alphaact);
Bfuzzy = vert2mu(Bg,alphaact);
dx  = Afuzzy*x+Bfuzzy*u;


%___________________________________
function fe = f_input(t,x)
% disturbio de entrada
% x = [x1 x2 t1 t2]
if t < 2 || t > 30
    fe = x;
else
    fe = x*1.2;
end


%_________________________________________________
%function [Af] = vert2mu(M)
function Af = vert2mu(M,alphas,iss)
%
vertices = length(M);
if nargin < 3
   iss = 1:vertices; 
end
%
A = zeros(size(M{1}));

j = 0;
for i = iss
    j = j + 1;
    A = A + alphas(1,j)*M{i};    
end
Af = A;


%_________________________________________________
function hiz = mshipfct(z,a,b)
% function hiz = mshipfct_ex2_FSS(z,a,b)
% Membership functions
% Inputs    -> z          : z = [z1 z2 z3 z4], 
%             a,b         : parameters
% Saidas   -> hiz : [h1 h2 h3 h4]
%
if nargin < 2
    a = 1.4;
    b = 0.7;
end
if z(1) > a
    z(1) = a;
    disp('x1 out of range!');
elseif z(1) < -a
    z(1) = -a; 
    disp('x1 out of range!');
end
if z(2) > b
    z(2) = b;
    disp('x3 out of range!');
elseif z(2) < -b
    z(2) = -b; 
    disp('x3 out of range!');
end

mu11 = z(1)^2/a^2;
mu12 = 1-mu11;
if z(2) == 0
    mu21 = 1;
else
    mu21 = (b*sin(z(2))-z(2)*sin(b))/(z(2)*(b-sin(b)));
end
mu22 = 1-mu21;
h1 = mu11*mu21; h2 = mu11*mu22; h3 = mu12*mu21; h4 = mu12*mu22;
hiz = [h1 h2 h3 h4];


%_____________________________________________________________
%function [P] = vert2mu_P(Pg)
function P = vert2mu_P(Pg,alphas,iss)
%
vertices = length(Pg);
if nargin < 3
   iss = 1:vertices; 
end
%
P = zeros(size(Pg{1}));
j = 0;
for i = iss
    j = j + 1;
    P = P + alphas(1,j)*Pg{i};
end

%_____________________________________________________________
%function [M] = Mx(M,x)
function M = Mx(Mg,x)
%
M = vert2mu_P(Mg,mshipfct_ex2_FSS(x));