function Ai = give_vertices
%give_vertices

[wf, V, Eref, m, n, Ro, Pref, Qref] = parameters;

global X
X = dxdt_fsolve;

z_lim = bounds_membership(X)

Ai = vertices(z_lim);