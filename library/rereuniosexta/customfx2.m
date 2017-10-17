function customfx2

a1=8;a2=20; b1=24; b2=25;
x1=-a1:0.1:a2;
x2=-b1:0.1:b2;
[X1,X2]=meshgrid(x1,x2);

% phase retrait nonlinear e.d.o
% dX = F(X)
F1=X2+sin(X1);
F2=16-(X1-5).^2;

figure(1);clf
streamslice(X1,X2,F1,F2,2);
xlabel('x_1');ylabel('x_2');title('Phase Retrait - nonlinear system');
line([1,9],[-sin(1),-sin(9)],'marker','o','linestyle','none','markerfacecolor','r')

x0 = [11; 10];
[t,y]=ode45(@fx,[0 50],x0); 
x02 = [16; 10];
[t2,y2]=ode45(@fx,[0 50],x02); 

hold on
plot(y(:,1),y(:,2),'r',y(1,1),y(1,2),'*r');
plot(y2(:,1),y2(:,2),'r',y2(1,1),y2(1,2),'*r');
hold off
axis([-a1 a2 -b1 b2])

%return

%Linearization
syms x1 x2
A = jacobian([x2+sin(x1); 16-(x1-5).^2],[x1 x2]);
x1=9; x2=-sin(9); %linearization point
%Ax = eval(subs(A));
Ax = subs(A);
eig(Ax)
F1x=Ax(1,1)*X1 + Ax(1,2)*X2;
F2x=Ax(2,1)*X1 + Ax(2,2)*X2;
figure(2);clf
streamslice(X1+x1,X2+x2,F1x,F2x,2);
xlabel('x_1');ylabel('x_2');title('Phase Retrait - linearized system');
line([1,9],[-sin(1),-sin(9)],'marker','o','linestyle','none','markerfacecolor','r')

return

function dx = fx(t,X)
F1=X(2)+sin(X(1));
F2=16-(X(1)-5).^2;
dx = [F1 F2]';
