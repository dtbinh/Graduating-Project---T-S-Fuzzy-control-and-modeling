%Curvas de Nivel
%V(x)=x'*P*x=gamma
close all; clc; clear all;

%ex=1 => Plota a curva de nível V(x)=x'*P*x=gamma, dados P e gamma
%ex=2 => Plota curvas 3-D de V(x)=x'*P*x=gamma, dados P e gamma
%ex=3 => Calcula os parametros de P para traj. conhecida, dados os ptos da trajetória
ex=1;

switch ex

    case 1
        %Plota a curva de nível V(x)=x'*P*x=gamma, dados P e gamma
        %a>0; c>0; a*c>b^2
        a=4;b=2;c=3; %a=2;b=2;c=3;
        gamma=1;

        P=[ a b;
            b c];
        dP=det(P);

        [V,D] = eig(P);
        v1=V(:,1);
        v2=V(:,2);
        %figure;plot(v1(1),v1(2),'r*',v2(1),v2(2),'b*');
        %hold on

        %Calculo
        x1max=sqrt(c*gamma/dP);
        x1=-x1max:1e-5:x1max;
        x2a=(-b*x1+sqrt(c*gamma-x1.^2*dP))/c;
        x2b=(-b*x1-sqrt(c*gamma-x1.^2*dP))/c;

        %Plot
        hold on
        plot(x1,x2a,'k')
        plot(x1,x2b,'k')
        hold off
        grid
        volume=det(inv(P)*gamma)
        
    case 2
        
        a=4;b=2;c=3; %a=2;b=2;c=3;
        gamma=1;
        P=[ a b;
            b c];
        mm=1;
        [X,Y] = meshgrid(-mm:2*mm/50:mm, -mm:02*mm/50:mm);
        Z = gamma*((X.^2).*a + (Y.^2).*c + (X.*Y).*b);
        MM=max(max(Z));
        [px,py] = gradient(Z,.2,.2);
        figure;contour(Z,50),hold on, quiver(px,py), hold off
        figure;[C,h] = contour(X,Y,Z);
            text_handle = clabel(C,h);
            set(text_handle,'BackgroundColor',[1 1 .6],...
            'Edgecolor',[.7 .7 .7])

        
        %figure;mesh(Z);
        %figure;mesh(Z);view(0,0);
        %figure;meshc(X,Y,Z);
        %axis([-mm mm -mm mm 0 MM])
        %figure;meshz(X,Y,Z)
        %axis([-mm mm -mm mm 0 MM])
        figure;surfc(X,Y,Z)
        axis([-mm mm -mm mm 0 MM])
        colormap hsv
        
    case 3
        %Calcula os parametros de P para traj. conhecida
        %(+-p1,0) ; (0,+-p2) ; (p3,p4)
        p1=0.22; p2=0.49; p3=0.164; p4=0.346;
        gamma=1;
        a=gamma/p1^2
        c=gamma/p2^2
        b=(gamma-a*p3^2-c*p4^2)/(2*p3*p4)

        % x1=[0.164; 0.346];
        % x2=[0.22; 0];
        % x3=[0; 0.49]
        % P=[ a b;
        %     b c];
        % x1'*P*x1

        P=[ a b;
            b c];
        dP=det(P);

        [V,D] = eig(P);
        v1=V(:,1);
        v2=V(:,2);
        %figure;plot(v1(1),v1(2),'r*',v2(1),v2(2),'b*');
        %hold on

        %Calculo
        x1max=sqrt(c*gamma/dP);
        x1=-x1max:1e-5:x1max;
        x2a=(-b*x1+sqrt(c*gamma-x1.^2*dP))/c;
        x2b=(-b*x1-sqrt(c*gamma-x1.^2*dP))/c;

        %Plot
        hold on
        plot(x1,x2a,'k')
        plot(x1,x2b,'k')
        hold off
        grid
        volume=det(inv(P)*gamma)
end



%Função para o calculo do volume do elipsoide
%exemplo: >>vol(2,2,3,1)
vol=@(a,b,c,gamma)det(inv([ a b; b c]))*gamma; %a>0; c>0; a*c>b^2; gamma>0