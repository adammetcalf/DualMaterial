close all;
clear;
clc;

syms L1 L2 L3 m1 m2 m3 Theta1(t) Theta2(t) Theta3(t) t g x1 y1 x2 y2 x3 y3;

x1 = L1*sin(Theta1);
x2 = x1 + L2*sin(Theta2);
x3 = x2 + L3*sin(Theta3);

y1 = -L1*cos(Theta1);
y2 = y1 - L2*cos(Theta2);
y3 = y2 - L3*cos(Theta3);

dotx1 = diff(x1,t);
doty1 = diff(y1,t);

dotx2 = diff(x2,t);
doty2 = diff(y2,t);

dotx3 = diff(x3,t);
doty3 = diff(y3,t);


v1Square = dotx1^2+doty1^2;
v1Square = expand(v1Square);
v1Square = simplify(v1Square);

v2Square = dotx2^2+doty2^2;
v2Square = expand(v2Square);
v2Square = simplify(v2Square);

v3Square = dotx3^2+doty3^2;
v3Square = expand(v3Square);
v3Square = simplify(v3Square);

% Form Kinetic Energies
sym T; 

T = (1/2)*m1*v1Square + (1/2)*m2*v2Square + (1/2)*m3*v3Square;
T = expand(T);
T = simplify(T);

% Form Potential Energies
syms V Vstiff k;

V = m1*g*y1 + m2*g*y2 + m3*g*y3;
V = expand(V);
V = simplify(V);

Vstiff = (1/2)*k*(Theta2-Theta1)^2;
Vstiff = expand(Vstiff);
Vstiff = simplify(Vstiff);

% Form Lagrangian
sym L;
L = T-(V);
L = expand(L);
L = simplify(L);
 
dL_wrt_dotT1 = diff(L,diff(Theta1));
dL_wrt_dotT1 = expand(dL_wrt_dotT1);
dL_wrt_dotT1 = simplify(dL_wrt_dotT1);

dL_wrt_dotT2 = diff(L,diff(Theta2));
dL_wrt_dotT2 = expand(dL_wrt_dotT2);
dL_wrt_dotT2 = simplify(dL_wrt_dotT2);

dL_wrt_dotT3 = diff(L,diff(Theta3));
dL_wrt_dotT3 = expand(dL_wrt_dotT3);
dL_wrt_dotT3 = simplify(dL_wrt_dotT3);

dt_wrt_dL_wrt_dotT1 = diff(dL_wrt_dotT1,t);
dt_wrt_dL_wrt_dotT1 = expand(dt_wrt_dL_wrt_dotT1);
dt_wrt_dL_wrt_dotT1 = simplify(dt_wrt_dL_wrt_dotT1);

dt_wrt_dL_wrt_dotT2 = diff(dL_wrt_dotT2,t);
dt_wrt_dL_wrt_dotT2 = expand(dt_wrt_dL_wrt_dotT2);
dt_wrt_dL_wrt_dotT2 = simplify(dt_wrt_dL_wrt_dotT2);

dt_wrt_dL_wrt_dotT3 = diff(dL_wrt_dotT3,t);
dt_wrt_dL_wrt_dotT3 = expand(dt_wrt_dL_wrt_dotT3);
dt_wrt_dL_wrt_dotT3 = simplify(dt_wrt_dL_wrt_dotT3);

dL_wrt_T1 = diff(L,Theta1);
dL_wrt_T1 = expand(dL_wrt_T1);
dL_wrt_T1 = simplify(dL_wrt_T1);

dL_wrt_T2 = diff(L,Theta2);
dL_wrt_T2 = expand(dL_wrt_T2);
dL_wrt_T2 = simplify(dL_wrt_T2);

dL_wrt_T3 = diff(L,Theta3);
dL_wrt_T3 = expand(dL_wrt_T3);
dL_wrt_T3 = simplify(dL_wrt_T3);

% calculate equations of motion
syms EQNmotionT1 EQNmotionT2

% Note, Eqns of motion are given by EQNmotionT1 = 0 and EQNmotionT2 = 0

EQNmotionT1 = dt_wrt_dL_wrt_dotT1 - dL_wrt_T1;
EQNmotionT1 = expand(EQNmotionT1);
EQNmotionT1 = simplify(EQNmotionT1);

EQNmotionT2 = dt_wrt_dL_wrt_dotT2 - dL_wrt_T2;
EQNmotionT2 = expand(EQNmotionT2);
EQNmotionT2 = simplify(EQNmotionT2);

EQNmotionT3 = dt_wrt_dL_wrt_dotT3 - dL_wrt_T3;
EQNmotionT3 = expand(EQNmotionT3);
EQNmotionT3 = simplify(EQNmotionT3);

syms a1 a2 a3 a4 a5 a6 f1 f2 f3 

A = [1 a1 a2; a3 1 a4; a5 a6 1];
F = [f1;f2;f3];
X = inv(A) * F