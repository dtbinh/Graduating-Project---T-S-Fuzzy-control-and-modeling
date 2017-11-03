function level_curve2(P, gamma, color)


x1_range =-pi/2:pi/2/1e3:pi/2;
x2_range = -pi/2:pi/2/1e3:pi/2;
[X1,X2] = meshgrid(x1_range,x2_range);
h1 = (1+X1)/2;
h2 = 1 - h1;
Pa11 = h1.*P{1}{1}(1,1)+ h2.*P{2}{1}(1,1);
Pa22 = h1.*P{1}{1}(2,2)+ h2.*P{2}{1}(2,2);
Pa12 = h1.*P{1}{1}(1,2)+ h2.*P{2}{1}(1,2);
Z = X1.*Pa11.*X1 + X2.*Pa22.*X2 + 2*X1.*Pa12.*X2;
%figure
%contour(X1,X2,Z,'ShowText','on')
v = [gamma, gamma];
hold on
% contour(X1,X2,Z,v,'LineColor',color);
contour(X1,X2,Z,'LineColor',color);
hold off

end