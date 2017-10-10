%clear; clc; close all

V=311;
Eref=311;
m=3.77/22000;
n=20/22000; %n=2000/22000 instavel
Ro=0.1;% se aumentar o valor de Ro fica mais lento e nao converge (delta=0.1 e n=20/22000).
wf=2*pi*60;
Pref=22000;
Qref=0;

Ts=1e-4;
Pf(1)=22000;
Qf(1)=0;
delta(1)=0.1;

for k=1:1000
    Pf(k+1)= (-wf*Pf(k)+(wf/Ro)*(V*(Eref-n*(Pf(k)-Pref))*cos(delta(k))-V^2))*Ts+Pf(k);
    Qf(k+1)=(-wf*Qf(k)-(wf*V*sin(delta(k))/Ro)*(Eref-n*(Pf(k)-Pref)))*Ts+Qf(k);
    delta(k+1)=(m*(Qf(k)-Qref))*Ts+delta(k);
end

figure
plot(Pf)
figure
plot(Qf)
figure
plot(delta)

