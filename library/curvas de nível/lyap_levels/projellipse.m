function vet = projellipse(P,linetype,couleur,holdflag)
% function projellisa(P,linetype,couleur,holdflag)
% visualisation de l'ellipsoide {x'*P*x <= 1} dans toutes les projections
% hold on;
if nargin == 1,
 linetype = '-';
 couleur = 'k';
end;
if nargin < 4,
 holdflag = 0;
end;

n = size(P,1);

if n == 2,
 nrow = 1; ncol = 1;
elseif n == 3,
 nrow = 3; ncol = 1;
elseif n == 4,
 nrow = 3; ncol = 2;
elseif n == 5,
 nrow = 5; ncol = 2;
elseif n == 6,
 nrow = 5; ncol = 3;
elseif n == 7,
 nrow = 7; ncol = 3;
else
 error('invalid dimension for P');
end;


x1 = 1; x2 = 2;
for i = 1:n*(n-1)/2,
 subplot(nrow, ncol, i);
 % -------------------------
 % Modifs pour la projection
 ind1 = [x1 x2];
 ind2 = [1:x1-1 x1+1:x2-1 x2+1:n];
 A = P(ind2, ind2);
 B = P(ind1, ind2);
 C = P(ind1, ind1);
 [h,vet]=plotellipse(C-B*inv(A)*B',zeros(n,1),linetype,couleur);
 % -------------------------
 xlabel(['x' int2str(x1)]);
 ylabel(['x' int2str(x2)]);
 if x2 < n,
  x2 = x2 + 1;
 else
  x1 = x1 + 1;
  x2 = x1 + 1;
 end;  
end;