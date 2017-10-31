

K = [-0.2 2];
G = 0.5*K;%[-0.3 4];
KG = K-G;
rho=1;

x2K1 = @(x1)(rho-K(1)*x1)/K(2);
x2K2 = @(x1)(-rho-K(1)*x1)/K(2);
x2KG1 = @(x1)(rho-KG(1)*x1)/KG(2);
x2KG2 = @(x1)(-rho-KG(1)*x1)/KG(2);

figure;
% Região com saturação
xl = 4;
fplot(x2K1,[-xl xl],'b')
hold on
fplot(x2K2,[-xl xl],'b')
% Região de setor
fplot(x2KG1,[-xl xl],'r')
hold on
fplot(x2KG2,[-xl xl],'r')
grid



% Calcul des matrices Xj constituees soit des element mu_i_sup soit des
% elements -mu_i_inf
% 2^n combinaisons de 1 ou -1
ns = 2;
binp = zeros(ns,2^ns);
binn = zeros(ns,2^ns);
for i = 1:ns
    for k = 0:2^i-1
        for j = k*(2^ns)/(2^i)+1:(k+1)*(2^ns)/(2^i)
            if (k/2 - fix(k/2)) == 0 binp(i,j) = 1; binn(i,j)=0;
            else                     binp(i,j) = 0; binn(i,j)=-1;
            end;
        end; % end for j
    end; % end for k
end;  %end for i

% On initialise les iterations avec le vecteur muX:
%muX = ones(2*n,1);
%muX = [5.0150;6.3978;7.0236;1.0934];
muX = [2.2 1.2 2.2 1.2]';

% Calcul des Xi
n=2;
for i = 1:2^n
    X(:,i) = [diag(binp(:,i)) diag(binn(:,i))]*muX;
end

% Plot box 2D
H=line([X(1,1);X(1,3);X(1,4);X(1,2);X(1,1)],[X(2,1);X(2,3);X(2,4);X(2,2);X(2,1)]);
set(H,'Color','magenta');

%Plot ellipse
gamma = 0.57;
P =[ 0.7120    -0.9621;
    -0.9621    2.6624];
gamma2 = 0.4;
P2 =[ 0.7120    0.9621;
    -0.9621    2.6624];
%projellisa(P*gamma,'k','-')    
vet{1} = projellipse(gamma*P);
vet{2} = projellipse(gamma2*P2);
[vet_inters,areapoly] = polyintersection(vet{1},vet{2},'c');
hold off;
box;
%gtext('U')
% print -depsc2 -r600 sat_regs
