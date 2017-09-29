x1=-4:.1:4;
x2=-5:.1:5;
[X1,X2]=meshgrid(x1,x2);

%A=diag([-1,-2])
%A=[-.5,-1;1,-.5]
A=inv([1,1;1,-1])*diag([1,-2])*[1,1;1,-1]

F1=A(1,1)*X1+A(1,2)*X2;
F2=A(2,1)*X1+A(2,2)*X2;

figure(1);clf
streamslice(X1,X2,F1,F2,.5);
axis([min(x1),max(x1),min(x2),max(x2)])
xlabel('x_1');ylabel('x_2');
line([0],[0],'marker','o','linestyle','none','markerfacecolor','r')
