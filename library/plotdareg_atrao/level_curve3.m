function level_curve3

P1 = [7 3;3 4];
%P2 = [1.6 1; 1 1.8];

th = pi/4;
R = [cos(th) -sin(th);sin(th) cos(th)];
P2 = R*P1;

P = {P1, P2, P2};

figure;
hold on
for i = 1:2
    x1 = -1:0.01:1;
    x2 = -1:0.01:1;
    [X1,X2] = meshgrid(x1,x2);
    Z = X1.*P{i}(1,1).*X1 + X2.*P{i}(2,2).*X2 + 2*X1.*P{i}(1,2).*X2;
    %figure
    %contour(X1,X2,Z,'ShowText','on')
    v = [1,1];
    contour(X1,X2,Z,v,'LineColor','b');
end
hold off

%generate the exponents of the monomials to create the variables P_k
vertices = 2;
deg = 2;
Kalpha = gen_coefs(vertices,deg);


%h1 = (1+cos(x1))/2;
%h2 = 1 - h1;
%figure;
hold on
x1 = -1:0.01:1;
x2 = -1:0.01:1;
[X1,X2] = meshgrid(x1,x2);
h1 = (1+X1)/2;
h2 = 1 - h1;
Pa11 = zeros(1); Pa22 = zeros(1); Pa12 = zeros(1);
for i = 1:length(Kalpha)
    Pa11 = Pa11 + h1.^Kalpha(i,2).*h2.^Kalpha(i,1).*P{i}(1,1);
    Pa22 = Pa22 + h1.^Kalpha(i,2).*h2.^Kalpha(i,1).*P{i}(2,2);
    Pa12 = Pa12 + h1.^Kalpha(i,2).*h2.^Kalpha(i,1).*P{i}(1,2);
end
Z = X1.*Pa11.*X1 + X2.*Pa22.*X2 + 2*X1.*Pa12.*X2;
%figure
%contour(X1,X2,Z,'ShowText','on')
v = [1,1];
contour(X1,X2,Z,v,'LineColor','r');
hold off

return
Pc{1}=P;
level_curve(Pc, 1, 'b')



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
