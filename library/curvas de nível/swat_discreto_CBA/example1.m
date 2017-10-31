function [Ai,Bi,vertices,order,inputs] = example1(b)

A1=[1   -b;
    -1  -0.5];
A2=[1   b;
    -1 -0.5];
Ai=[A1 A2];

B1=[5+b;
    2*b];
B2=[5-b;
    -2*b];
Bi=[B1 B2];

[order,vertices] = size(Ai);
vertices = vertices / order;
[order,inputs] = size(Bi);
inputs = inputs / vertices;