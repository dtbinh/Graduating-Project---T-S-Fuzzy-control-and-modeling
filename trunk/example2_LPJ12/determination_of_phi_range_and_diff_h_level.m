function determination_of_phi_range_and_diff_h_level(n, h, A, xi, x1, x2, phi)

A_nonlinear = 0;
for i = 1:n
    A_nonlinear = A_nonlinear + h(i)*A(:, :, i);
end

syms x1 x2
for i = 1:n
    dh(i) = diff(h(i))*(A_nonlinear(1, :) * [x1; x2]);
end


xx1 = xi(1):0.1:xi(2);
xx2 = xi(1):0.1:xi(2);
[X1,X2] = meshgrid(xx1,xx2);

for k = 1:n
    for i=1:length(xx1)
        for j=1:length(xx2)
            x1 = X1(i,j);
            x2 = X2(i,j);
            Zi(i,j) = subs(dh(k));
        end
    end
    Z{k} = Zi;
end

v = [phi,phi];
%figure;
hold on
contour(X1,X2,Z{1},v,'LineColor','b');
contour(X1,X2,Z{2},v,'LineColor','r');
hold off

end

