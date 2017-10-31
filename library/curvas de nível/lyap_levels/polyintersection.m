function [vet_out,areapoly] = polyintersection(v1,v2,color)
%function vet_out = polyintersection(v1,v2,fillpoly)
% Retorna o vetor da interseçao de 2 politopos centrados na origem
% necessario ter boa resolução dos vetores (> 300 elementos) para evitar distorções
% otimizado para 1000 elementos
% Entradas -> v1         : politopo1 (x,y) centrado na origem
%             v2         : politopo1 (x,y) centrado na origem
%             color      : cor do preenchimento, parametro opcional (ex. color='c')   
% Saidas   -> vet_out    : vetores da interseçao dos 2 politopos
%             areapoly   : area da região da interseçao dos 2 politopos
% Created: 06/mar/2010
% Author: edutog@dt.fee.unicamp.br
% Exemplo
% figure;
% P{1}=[0.0171 0.000783;0.000783 0.00278];
% P{2}=[0.00942 0.00392;0.00392 0.00558];
% for i=1:2
%     vet{i} = projellisa2b(P{i});
% end
% plot(vet{1}(:,1),vet{1}(:,2),'b');
% hold on
% plot(vet{2}(:,1),vet{2}(:,2),'r');
% [vet_inters,areapoly] = polyintersection(vet{1},vet{2},'c');
% hold off;box;grid;

[th1,rho1] = cart2pol(v1(:,1),v1(:,2));
v1polar = [th1 rho1];
[th2,rho2] = cart2pol(v2(:,1),v2(:,2));
v2polar = [th2 rho2];

for i = 1:length(th1)-1
    dif(i) = th1(i+1)-th1(i);
end
diff1 = max(dif(i));
for i = 1:length(th2)-1
    dif(i) = th2(i+1)-th2(i);
end
diff2 = max(dif(i));
delta=max(diff1,diff2);
%delta = 0.01; %0.2 (se o no. de elementos dos vetores for pequeno aumentar o valor de delta)

vet_inters = [];
for th = -pi:delta:pi+delta
   [rows1] = find(v1polar(:,1)>=th & v1polar(:,1)<=th+delta);
   [rows2] = find(v2polar(:,1)>=th & v2polar(:,1)<=th+delta);
   if isempty(rows1) || isempty(rows2)
      %disp('empty');
   else
   if max(v1polar(rows1,2)) > max(v2polar(rows2,2))
       vet_inters = [vet_inters; v2polar(rows2,:)];
   else       
       vet_inters = [vet_inters; v1polar(rows1,:)];
   end
   end
end
[x,y] = pol2cart(vet_inters(:,1),vet_inters(:,2));
vet = [x y]; vet_out = vet;
if nargin >= 3
    fill(vet(:,1),vet(:,2),color);
    hold on
end
plot(vet(:,1),vet(:,2),'k');
%hold off;
box;%grid;
areapoly=polyarea(vet(:,1),vet(:,2));
