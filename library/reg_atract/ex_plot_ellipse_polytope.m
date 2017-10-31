clear all;

n=2;
muX = [1.5;1.5];%ones(n,1); 

% Calcul des matrices Xj constituees soit des element mu_i_sup soit des
% elements -mu_i_inf
% 2^n combinaisons de 1 ou -1
ns = n;
binpn = zeros(ns,2^ns);
for i = 1:ns
    for k = 0:2^i-1
        for j = k*(2^ns)/(2^i)+1:(k+1)*(2^ns)/(2^i)
            if (k/2 - fix(k/2)) == 0 binpn(i,j) = 1;
            else                     binpn(i,j) = -1;
            end;
        end; % end for j
    end; % end for k
end;  %end for i

% Calcul des Xi
for i = 1:2^n
    X(:,i) = diag(binpn(:,i))*muX;
    % diag(v) returns a square diagonal matrix with the elements of vector
    % v on the main diagonal.
end

% plots x_k
line([X(1,1);X(1,3);X(1,4);X(1,2);X(1,1)],[X(2,1);X(2,3);X(2,4);X(2,2);X(2,1)]);

%return

A=rand(2,2);
A = A - eye(2)*max([0 max(eig(A))+0.5]);
eig(A)
In = eye(n);
m = 1;

opt = 'P';
% opt = 'W';

if opt == 'P'

    P = sdpvar(n,n,'symmetric');
    gamma = sdpvar(1,1,'symmetric');
    T = sdpvar(n,n,'symmetric');
    
    LMIs = [];
    LMIs = [LMIs, P >= 0];
    LMIs = [LMIs, A'*P + P*A <= 0];
    
    for i = 1:n
        MatQ = [P*muX(i) In(i,:)'; In(i,:) gamma*muX(i)];
        LMIs = [LMIs, MatQ >= 0];
    end
    LMIs = [LMIs, [T P;P P]>=0];
    
    bet1 = 1;
    bet2 = 1;
    crit = bet1*trace(T)+bet2*gamma;
    %crit = [];
    
    solvesdp(LMIs,crit,sdpsettings('solver','sedumi','verbose',0));
    [p,d]=checkset(LMIs);
    pmin = min(checkset(LMIs));
    display(pmin)
    
    tol = 1e-7;
    if sum(p < -tol)
        msgbox 'instavel'
    end
    P=double(P);
    W = inv(P);
    gamma=double(gamma);
    
    projellisa(P*gamma,'b','-')
    
    volume = sqrt(det(W/gamma))
    
    
else
    W = sdpvar(n,n,'symmetric');
    gamma = sdpvar(1,1,'symmetric');
    T = sdpvar(n,n,'symmetric');
    
    LMIs = [];
    LMIs = [LMIs, W >= 0];
    LMIs = [LMIs, W*A' + A*W <= 0];
    
    for i = 1:n
        MatQ = [muX(i)*W W*In(i,:)';In(i,:)*W gamma*muX(i)];
        LMIs = [LMIs, MatQ >= 0];
    end
    LMIs = [LMIs, [T In;In W]>=0];
    
    bet1 = 1;
    bet2 = 1;
    crit = bet1*trace(T)+bet2*gamma;
    %crit = [];
    
    solvesdp(LMIs,crit,sdpsettings('solver','sedumi','verbose',0));
    [p,d]=checkset(LMIs);
    pmin = min(checkset(LMIs));
    display(pmin)
    
    tol = 1e-7;
    if sum(p < -tol)
        msgbox 'instavel'
    end
    W=double(W);
    P = inv(W);
    gamma=double(gamma);
    
    projellisa(P*gamma,'b','-')
    
    volume = sqrt(det(W/gamma))
    
end