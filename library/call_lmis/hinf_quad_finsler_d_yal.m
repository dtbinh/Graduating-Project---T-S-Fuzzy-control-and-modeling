function output = hinf_quad_finsler_d_yal(Ai,B1i,B2i,Ci,D1i,D2i,param)
% function output = pid_hinf_quad_finsler_d_yal(Ai,B1i,B2i,Ci,D1i,D2i,param)
%
% Synthesize a robust state feedback control gain that ensure closed-loop quadratic
% stability with a prescibed H-infinity attenuation level for discrete-time time-varying polytopic system.
% The LMIs are programmed using YALMIP  and can be solved  by any LMI solver supported by
% YALMIP (SeDuMi is the default).
%
% input: Ai,B1i,B2i,Ci,D1i,D2i -> vertices of the system matrices.
%        param         -> optional parameters
%        param.hinf    -> if this parameter is different of zero, the routine will
%                                 check the feasibility of the LMIs for this value of the hinf norm.
%        param.maxViolation -> maximum violation of the constraints admissible by the
%                              program to consider the obtained guaranteed cost still valid. 1e-7 is the
%                              default.
%
% output:output.hinf     -> h-infinity guaranteed cost (0 if unfeasible)
%        output.cpusec   -> cpu time to solve the LMIs (seconds)
%        output.cpusec_m -> cpu time to mount the LMIs (seconds)
%        output.K        -> number of scalar variables
%        output.L        -> number of LMIs rows
%        output.W        -> solution variable W (Lyapunov matrix)
%        output.Ki       -> synthesized gains ([Z_1 Z_2 .. Z_N]*inv(W)]
%        output.delta    -> minimal primal residual returned by the LMI solver (SeDuMi is the default).
%
% Example: 2 states and 2 vertices
%
% A1 = [-0.9  0.2;
%       -0.5  -1.9 ];
% A2 = [-1.1  -0.1;
%       0.5    -2.8];
% B1 = [1;0];  B2 = [0;-1];
% Ai = [A1 A2]; B1i = [B1 B2];
% B2i = B1i; Ci = [1 0 1 0;0 1 0 1];
% D1i = zeros(2,2); D2i = [0.3 -0.2;0.2 -0.2];
% Info = ssf_hinf_quad_d_yal(Ai,B1i,B2i,Ci,D1i,D2i)
%
% Date: 01/08/2016
% Author: estognetti@ene.unb.br

%determine the number of vertices of polytope
[order,vertices] = size(Ai);
vertices = vertices / order;
inputs_c = size(B2i,2)/vertices;
inputs_d = size(B1i,2)/vertices;
outputs = size(Ci,1);

if nargin == 7
    if isfield(param,'struct')
        struct = param.struct;
    else
        struct = 0;
    end
    if isfield(param,'hinf')
        hinf = param.hinf;
    else
        hinf = 0;
    end
    if isfield(param,'max_gain')
        max_gain = param.max_gain;
    else
        max_gain = 0;
    end
    if isfield(param,'maxViolation')
        maxViolation = param.maxViolation;
    else
        maxViolation = 1e-7;
    end

else
    struct = 0;
    hinf = 0;
    max_gain = 0;
    maxViolation = 1e-7;
end

output.cpusec_m = clock;

%new LMI system
LMIs = [];

if hinf == 0
    sigma = sdpvar(1);
    obj = sigma;
else
    sigma = hinf*hinf;
    obj = [];
end

switch struct
    case 0 % full structure
         W = sdpvar(order,order,'symmetric');
         Z = sdpvar(inputs_c,order,'full');            
         G = sdpvar(order,order,'full');
        
    case 1 % blocked structure
        W = sdpvar(order,order,'symmetric');
        g = []; z = [];
        for i = 1:order
            g = [g sdpvar(1)];
            if inputs_c == order
                z = [z sdpvar(1)];
            end
        end
        G = diag(g);
        if ~isempty(z), Z = diag(z);else Z = sdpvar(inputs_c,order,'full'); end
        
end

% LMIs
LMIs = [LMIs, W > 0];
for i=1:vertices
    A  = Ai(:, (i-1)*order+1   :i*order);
    C  = Ci(:, (i-1)*order+1   :i*order);
    B1 = B1i(:,(i-1)*inputs_d+1:i*inputs_d);
    D1 = D1i(:,(i-1)*inputs_d+1:i*inputs_d);
    B2 = B2i(:,(i-1)*inputs_c+1:i*inputs_c);
    D2 = D2i(:,(i-1)*inputs_c+1:i*inputs_c);


    T11 = W;
    T12 = A*G + B2*Z;
    T13 = zeros(order,outputs);
    T14 = B1;
    T22 = G+G'-W;
    T23 = G*C' +Z'*D2';
    T24 = zeros(order,inputs_d);
    %T33 = sigma*eye(outputs);
    T33 = eye(outputs);
    T34 = zeros(outputs,inputs_d);
    %T44 = eye(inputs_d);
    T44 = sigma*eye(inputs_d);

    % Stability + hinf
    LMIs = [LMIs,     [T11  T12  T13  T14;
                       T12' T22  T23  T24;
                       T13' T23' T33  T34;
                       T14' T24' T34' T44] > 0];

% Stability only
% MIs = [LMIs, [T11  T12; T12' T22] > 0];
                   
end

if max_gain > 0
    %MIs = [LMIs, [W Z';Z max_gain*max_gain*eye(inputs_c)] > 0];
end

% determine the number of LMI rows
output.L = 0;
for i=1:size(LMIs,1)
    output.L = output.L + size(LMIs{i},1);
end
% determine the number of scalar variables
output.V = size(getvariables(LMIs),2);

% evaluate the elapsed time to mount the LMIs set
output.cpusec_m = etime(clock,output.cpusec_m);

% solve the LMIs
sol = solvesdp(LMIs,obj,sdpsettings('verbose',0,'solver','sedumi'));
%sol = solvesdp(LMIs,obj,sdpsettings('verbose',0,'solver','lmilab'));
%sol = solvesdp(LMIs,obj,sdpsettings('verbose',0,'solver','mosek'));

% evaluate the elapsed time to solve the LMIs set
output.cpusec = sol.solvertime;

% retrieving the minimal primal residual
p=min(checkset(LMIs));
output.delta = p;

output.hinf = 0;
%capturing the solutions (if ones exist)
if p  > -maxViolation %adopted precision for the minimum primal residual
    output.W = double(W);
    output.G = double(G);
    output.Z = double(Z);
    output.K = double(Z)*inv(double(G));
    output.hinf = sqrt(double(sigma));
end

