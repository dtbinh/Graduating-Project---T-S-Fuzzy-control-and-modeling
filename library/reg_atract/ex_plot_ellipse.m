
P=rand(2,2);
P = P - eye(2)*min([0 min(eig(P))-0.1]);
P=P+P';
eig(P)
gamma=1;

% function projellisa(P,linetype,couleur,holdflag)
% visualisation de l'ellipsoide {x'*P*x <= 1} dans toutes les projections
projellisa(P*gamma,'b','-')

% function h=plotellip(A,b,linetype,couleur)
% [xe,ye]=plotellip(A,xc,linetype);
%   plots ellipsoid (x-xc)*A*(x-xc) = 1 in R^2, or
% [xe,ye]=plotellip(A,b,linetype);
%   plots ellipsoid x'*A*x + b'*x + c = 0 in R^2
plotellisa(P,rand(2,1),'k','-')

