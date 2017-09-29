x1=-5:.1:5;
x2=-5:.1:5;
[X1,X2]=meshgrid(x1,x2);

% asympto. stable
%A1=[-1,.25;-1,-1]
%A2=[-1,-1;.25,-1]

% asympto. stable but no common quadratic Lyapunov function
A1=[-1,-1;1,-1]
A2=[-1,-8;1/8,-1]

% just stable
%A1=[0,-1;1,-2]
%A2=[0,1;-1,-2]

eig(A1+A1')
eig(A2+A2')


F1=A1(1,1)*X1+A1(1,2)*X2;
F2=A1(2,1)*X1+A1(2,2)*X2;



figure(1);clf
streamslice(X1,X2,F1,F2,.25);
axis([min(x1),max(x1),min(x2),max(x2)])
axis('square')
xlabel('x_1');ylabel('x_2');
line([0],[0],'marker','o','linestyle','none','markerfacecolor','r')
