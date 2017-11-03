function level_curve_v(P_n, gamma, xi, h, c)

syms x1 x2
n=2;
%generate the exponents of the monomials to create the variables P_k
vertices = 2;
deg = 1;
Kalpha = gen_coefs(vertices,deg);

xx1 = xi(1):0.1:xi(2);
xx2 = xi(1):0.1:xi(2);
[X1,X2] = meshgrid(xx1,xx2);

for k = 1:n
    for i=1:length(xx1)
        for j=1:length(xx2)
            x1 = X1(i,j);
            x2 = X2(i,j);
            hh(i,j) = subs(h(k));
        end
    end
    H{k} = hh;
end

Pa11 = zeros(1); Pa22 = zeros(1); Pa12 = zeros(1);
for i = 1:length(Kalpha)
    Pa11 = Pa11 + H{1}.^Kalpha(i,2).*H{2}.^Kalpha(i,1).*P_n{i}{1}(1,1);
    Pa22 = Pa22 + H{1}.^Kalpha(i,2).*H{2}.^Kalpha(i,1).*P_n{i}{1}(2,2);
    Pa12 = Pa12 + H{1}.^Kalpha(i,2).*H{2}.^Kalpha(i,1).*P_n{i}{1}(1,2);
end
Z = X1.*Pa11.*X1 + X2.*Pa22.*X2 + 2*X1.*Pa12.*X2;

v = [gamma,gamma];
hold on
contour(X1,X2,Z,v,'LineColor',c);
hold off



%___________________________________________________________________________________________________
function [comb] = gen_coefs(vertices,degree)
%generates all the solutions of: a_1 + a_2 + ... + a_vertices = degree

name = strcat('comb_N',int2str(vertices),'d',int2str(degree),'.mat');
if exist(name,'file') == 2
    load(name);
else
    comb = [];
    if vertices == 1
        comb = degree;
    else
        for i = 0:degree
            temp = gen_coefs(vertices-1,degree-i);
            comb   = [ comb
                i*ones(size(temp,1),1) temp ];
        end
    end
end
