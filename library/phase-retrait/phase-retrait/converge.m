x1=-2:.01:2;
x2=-3:.01:3;
[X1,X2]=meshgrid(x1,x2);

F1=X1.^2-X2.^2;
F2=2.*X1.*X2;

hold off
figure(1);clf
streamslice(X1,X2,F1,F2,2);
xlabel('x_1');ylabel('x_2');
line([0],[0],'marker','o','linestyle','none','markerfacecolor','r')
hold off