x1=-5:.1:5;
x2=-5:.1:5;
[X1,X2]=meshgrid(x1,x2);


A1=[-.5,-.4;3,-.5]
A2=[-.5,-3;.4,-.5]


F1=A1(1,1)*X1+A1(1,2)*X2;
F2=A1(2,1)*X1+A1(2,2)*X2;

k=find(X1.*X2<=0);

F1(k)=A2(1,1)*X1(k)+A2(1,2)*X2(k);
F2(k)=A2(2,1)*X1(k)+A2(2,2)*X2(k);


figure(1);clf
streamslice(X1,X2,F1,F2,.25);
axis([min(x1),max(x1),min(x2),max(x2)])
axis('square')
xlabel('x_1');ylabel('x_2');
line([0],[0],'marker','o','linestyle','none','markerfacecolor','r')
