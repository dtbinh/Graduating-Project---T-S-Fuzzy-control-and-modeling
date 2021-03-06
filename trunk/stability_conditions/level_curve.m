function level_curve(P)
    %Plots the level curve V(x) = x'*P*x = gamma, given P and gamma.
        
    gamma=1;

    for i = 1:length(P)
        Pi = P{i}{1};
        dP = det(Pi);

        [V,D] = eig(Pi);
        %v1=V(:,1);
        %v2=V(:,2);
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
        %volume = det(inv(Pi)*gamma) 
    end
end

