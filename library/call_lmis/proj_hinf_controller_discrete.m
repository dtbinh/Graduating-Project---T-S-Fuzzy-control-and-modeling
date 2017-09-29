clear all; close all; clc;
% author: Eduardo S Tognetti <estognetti@ene.unb.br>
% date: 03/09/2016

imprimir = 0;
if imprimir == 1, diary diaryout.txt; format long; else format short; end

% Dados
load dados1 % Sistema discreto
%load dados2 % Sistema discreto c/ 2 vértices na forma [A1 A2], ...
Ts = 1;

% LMI-based hinf state feedback controller
param.struct = 0; %0 (full) or 1 (blocked)
%param.hinf = 6; norma Hinf dada

out = hinf_quad_finsler_d_yal(Ai,Bid,Bi,Ci,Did,Di,param)
    

%determine the number of vertices of polytope
[order,vertices] = size(Ai);
vertices = vertices / order;
inputs_c = size(Bi,2)/vertices;
inputs_d = size(Bid,2)/vertices;
outputs = size(Ci,1);
for i=1:vertices
    A{i}  = Ai(:, (i-1)*order+1   :i*order);        
    Bw{i} = Bid(:,(i-1)*inputs_d+1:i*inputs_d);
    Bu{i} = Bi(:,(i-1)*inputs_c+1:i*inputs_c);    
    C{i}  = Ci(:, (i-1)*order+1   :i*order);
    Dw{i} = Did(:,(i-1)*inputs_d+1:i*inputs_d);
    Du{i} = Di(:,(i-1)*inputs_c+1:i*inputs_c);
    K = out.K; %robust gain
    
    % Malha fechada sist. discreto com controlador projetado
    Amfd{i} = A{i} + Bu{i}*K;
    display(eig(Amfd{i}))

    % Plot: MF controlador projetado
    sysmfdd = ss(Amfd{i},Bw{i},C{i}+Du{i}*K,Dw{i},Ts);
    figure; step(sysmfdd);
    if imprimir == 1, print -depsc steptraject; end
    
    figure; bode(sysmfdd)
    if imprimir == 1, print -depsc bodeaug; end
end

if imprimir == 1, diary OFF; end
 

%return

