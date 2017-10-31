function Info = synth_swsat_tais_slack_poleplace_otim_d_yal(Ai,Bi,rhoi,c,r)
%function Info = synth_swsat_tais_slack_poleplace_otim_d_yal(Ai,Bi,rhoi,c,r)
% inputs: Ai            -> Ai = [A1 A2 ... AN]
%         Bi            -> Bi = [B1 B2 ... BN]
%         rhoi          -> u(t) = rhoi, |Kx| > rhoi 
%
% outputs: Info.feas     -> stable (1) or not (0)
%          Info.cpusec   -> cpu time to solve the LMIs (seconds)
%          Info.cpusec_m -> cpu time to mount the LMIs (seconds)
%          Info.Wi       -> Solution variables Wi
%          Info.Y        -> Solution variables Yi
%          Info.X        -> Solution variables X
%          Info.Zi       -> Solution variables Zi
%          Info.S        -> Solution variables S
%          Info.Ki       -> Controller Ki = Zi * inv(X)
%          Info.beta    -> menor eixo da regiao elipsoidal de estabilidade robusta do sistema chaveado
%          Info.K        -> number of scalar variables
%          Info.L        -> number of LMI rows
% 
% Example:
% b=1;
% A1=[1   -b;
%     -1  -0.5];
% A2=[1   b;
%     -1 -0.5];
% Ai=[A1 A2];
% 
% B1=[5+b;
%     2*b];
% B2=[5-b;
%     -2*b];
% Bi=[B1 B2];
% rhoi = [10];
% c = 0.2;
% r = 0.5;
% Info = synth_swsat_tais_slack_poleplace_otim_d_yal(Ai,Bi,rhoi,c,r)
%
% 08/fev/2010
% taisrc@dt.fee.unicamp.br
%
%determine the number of vertices of polytope
[order,vertices] = size(Ai);
vertices = vertices / order;
[order,inputs] = size(Bi);
inputs = inputs / vertices;

Info.cpusec_m = clock;

%new LMI system
LMIs = set([]);

Info.L = 0;

%Variables
for i = 1:vertices
    W{i} = sdpvar(order,order,'symmetric'); %W = W' > 0        
    Z{i} = sdpvar(inputs,order,'full');
    Y{i} = sdpvar(inputs,order,'full');
    X{i} = sdpvar(order,order,'full');
end

for i = 1:vertices
    sii{i} = [];
    for j = 1:inputs
        s(i,j) = sdpvar(1,1); %S = diag(sii) > 0
        sii{i} = [sii{i} s(i,j)];
    end
    S{i} = diag(sii{i});
end
beta = sdpvar(1,1);
Im=eye(inputs);
In=eye(order);

% LMIs - positive definite matrices
for i=1:vertices    
    %LMIs = LMIs + set(W{i} > 0);
    %LMIs = set([]);
    LMIs = [LMIs, W{i} >= 0];
    Info.L = Info.L + order;

    %LMIs = LMIs + set(S{i} > 0);
    LMIs = [LMIs, S{i} >= 0];
    Info.L = Info.L + order;
end

% LMIs - conditions for local stability
for j=1:vertices
    for i = 1:vertices

        Aii = Ai(:,order*(i-1)+1:i*order);
        Bii = Bi(:,inputs*(i-1)+1:i*inputs);
        
        T11 = r*(W{i} - X{i} - X{i}');
        T12 = X{i}'*Aii' + Z{i}'*Bii' + c*X{i}';
        T13 = Y{i}';
        T22 = -r*W{j};
        T23 = -Bii*S{i};
        T33 = -2*S{i};
        
        TT =[T11  T12  T13;
             T12' T22  T23;
             T13' T23' T33];
        %LMIs = LMIs + set(TT < 0);
        LMIs = [LMIs, TT <= 0];
        Info.L = Info.L + 2*order + inputs;

    end
end

% LMIs - conditions for ellipsoid belonging to polyhedral region
% (condicao para que elipsoide esteja contido na regiao poliedral)
for i=1:vertices
    for j=1:inputs
        T11 = W{i};
        T12 = Z{i}(j,:)' - Y{i}(j,:)';
        T22 = rhoi(j)^2;
        TT = [ T11 T12; T12' T22];        
        %LMIs = LMIs + set(TT >= 0);
        LMIs = [LMIs, TT >= 0];
        Info.L = Info.L + order + inputs;
    end 
end

% LMIs - conditions for minimization of maximum eigenvalue of W
% (condicao para minimizacao do maximo autovalor de W)
for i=1:vertices
        T11 = beta*In;
        T12 = In;
        T22 = W{i};
        TT = [ T11 T12; T12' T22];
        %LMIs = LMIs + set(TT >= 0);
        LMIs = [LMIs, TT >= 0];
        Info.L = Info.L + 2*order;    
end


%LMIs solutions - otimization problem
Info.cpusec_m = etime(clock,Info.cpusec_m);

Info.K = size(getvariables(LMIs),2);

sol = solvesdp(LMIs,beta,sdpsettings('verbose',0,'solver','sedumi'));
Info.cpusec = sol.solvertime;
[p,d]=checkset(LMIs);

Info.feas = 0;
%capturing the solutions (if ones exist)
if sum(p < 0) == 0
    Info.Si = [];
    Info.Yi = [];
    Info.Xi = [];
    Info.Wi = [];
    Info.Zi = [];
    Info.Ki = [];
    for i=1:vertices
        Info.Si = [Info.Si double(S{i})];        
        Info.Wi = [Info.Wi double(W{i})];
        Info.Zi = [Info.Zi double(Z{i})];
        Info.Yi = [Info.Yi double(Y{i})];
        Info.Xi = [Info.Xi double(X{i})];
        Info.Ki = [Info.Ki double(Z{i})*inv(double(X{i}))];
    end
    Info.beta = double(beta);
    Info.feas = 1;
end
