function phase_retrait(x)

plot3(x(:,1),x(:,2),x(:,3));
title ('Phase portrait' );
xlabel('Pf'); ylabel ('Qf' ); zlabel ('delta' );
line([x(end,1)],[0], [0],'marker','o','linestyle','none','markerfacecolor','r')
grid;
hold off

end

