x1=-4:.01:4;
x2=-5:.01:5;
[X1,X2]=meshgrid(x1,x2);

F1=X2;
F2=-X1-.5*(1-X1.^2).*X2;

figure(1);clf
streamslice(X1,X2,F1,F2,1);
axis([min(x1),max(x1),min(x2),max(x2)])
xlabel('x_1');ylabel('x_2');
line([0],[0],'marker','o','linestyle','none','markerfacecolor','r')
