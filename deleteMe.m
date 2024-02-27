close all;
clear;
clc;

syms l1 l2 l3 l4 m1 m2 m3 m4 T1 T2 T3 T4 T1t T2t T3t T4t

dotx1 = l1*T1t*cos(T1);
doty1 = l1*T1t*sin(T1);

dotx2 = dotx1 + l2*T2t*cos(T2);
doty2 = doty1 + l2*T2t*sin(T2);

dotx3 = dotx1 + dotx2 + l3*T3t*cos(T3);
doty3 = doty1 + doty2 + l3*T3t*sin(T3);

dotx4 = dotx1 + dotx2 + dotx3 + l4*T4t*cos(T4);
doty4 = doty1 + doty2 + doty3 + l4*T4t*sin(T4);

v1Square = dotx1^2+doty1^2;
v1Square = expand(v1Square);
v1Square = simplify(v1Square)

v2Square = dotx2^2+doty2^2;
v2Square = expand(v2Square);
v2Square = simplify(v2Square)

v3Square = dotx3^2+doty3^2;
v3Square = expand(v3Square);
v3Square = simplify(v3Square);

v4Square = dotx4^2+doty4^2;
v4Square = expand(v4Square);
v4Square = simplify(v4Square)
