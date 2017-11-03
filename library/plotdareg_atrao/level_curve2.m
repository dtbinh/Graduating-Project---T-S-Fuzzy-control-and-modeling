clear
P1 = [7 3;3 4];
%P2 = [1.6 1; 1 1.8];

th = pi/4;
R = [cos(th) -sin(th);sin(th) cos(th)];
P2 = R*P1;

P = {P1, P2};

figure;
hold on
for i = 1:2
x1 = -1:0.01:1;
x2 = -1:0.01:1;
[X1,X2] = meshgrid(x1,x2);
Z = X1.*P{i}(1,1).*X1 + X2.*P{i}(2,2).*X2 + 2*X1.*P{i}(1,2).*X2;
%figure
%contour(X1,X2,Z,'ShowText','on')
v = [1,1];
contour(X1,X2,Z,v)
end
hold off


%h1 = (1+cos(x1))/2;
%h2 = 1 - h1;
%figure;
hold on
x1 = -1:0.01:1;
x2 = -1:0.01:1;
[X1,X2] = meshgrid(x1,x2);
h1 = (1+X1)/2;
h2 = 1 - h1;
Pa11 = h1.*P1(1,1)+ h2.*P2(1,1);
Pa22 = h1.*P1(2,2)+ h2.*P2(2,2);
Pa12 = h1.*P1(1,2)+ h2.*P2(1,2);
Z = X1.*Pa11.*X1 + X2.*Pa22.*X2 + 2*X1.*Pa12.*X2;
%figure
%contour(X1,X2,Z,'ShowText','on')
v = [1,1];
contour(X1,X2,Z,v,'LineColor','r');
hold off

return
Pc{1}=P;
level_curve(Pc, 1, 'b')