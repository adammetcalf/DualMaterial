close all;
clear;
clc;

syms L1 L2 m1 m2 Theta1(t) Theta2(t) t g x1 y1 x2 y2;

x1 = L1*sin(Theta1);
x2 = x1 + L2*sin(Theta2);

y1 = -L1*cos(Theta1);
y2 = y1 - L2*cos(Theta2);

dotx1 = diff(x1,t);
doty1 = diff(y1,t);

dotx2 = diff(x2,t);
doty2 = diff(y2,t);


v1Square = dotx1^2+doty1^2;
v1Square = expand(v1Square);
v1Square = simplify(v1Square);

v2Square = dotx2^2+doty2^2;
v2Square = expand(v2Square);
v2Square = simplify(v2Square);

% Form Kinetic Energies
sym T; 

T = (1/2)*m1*v1Square + (1/2)*m2*v2Square;
T = expand(T);
T = simplify(T);

% Form Potential Energies
syms V Vstiff k;

V = m1*g*y1 + m2*g*y2;
V = expand(V);
V = simplify(V);

Vstiff = (1/2)*k*(Theta2-Theta1)^2;
Vstiff = expand(Vstiff);
Vstiff = simplify(Vstiff);

% Form Lagrangian
sym L;
L = T-(V+Vstiff);
L = expand(L);
L = simplify(L);
 
dL_wrt_dotT1 = diff(L,diff(Theta1));
dL_wrt_dotT1 = expand(dL_wrt_dotT1);
dL_wrt_dotT1 = simplify(dL_wrt_dotT1);

dL_wrt_dotT2 = diff(L,diff(Theta2));
dL_wrt_dotT2 = expand(dL_wrt_dotT2);
dL_wrt_dotT2 = simplify(dL_wrt_dotT2);

dt_wrt_dL_wrt_dotT1 = diff(dL_wrt_dotT1,t);
dt_wrt_dL_wrt_dotT1 = expand(dt_wrt_dL_wrt_dotT1);
dt_wrt_dL_wrt_dotT1 = simplify(dt_wrt_dL_wrt_dotT1);

dt_wrt_dL_wrt_dotT2 = diff(dL_wrt_dotT2,t);
dt_wrt_dL_wrt_dotT2 = expand(dt_wrt_dL_wrt_dotT2);
dt_wrt_dL_wrt_dotT2 = simplify(dt_wrt_dL_wrt_dotT2);

dL_wrt_T1 = diff(L,Theta1);
dL_wrt_T1 = expand(dL_wrt_T1);
dL_wrt_T1 = simplify(dL_wrt_T1);

dL_wrt_T2 = diff(L,Theta2);
dL_wrt_T2 = expand(dL_wrt_T2);
dL_wrt_T2 = simplify(dL_wrt_T2);

% calculate equations of motion
syms EQNmotionT1 EQNmotionT2

% Note, Eqns of motion are given by EQNmotionT1 = 0 and EQNmotionT2 = 0

EQNmotionT1 = dt_wrt_dL_wrt_dotT1 - dL_wrt_T1;
EQNmotionT1 = expand(EQNmotionT1);
EQNmotionT1 = simplify(EQNmotionT1)

EQNmotionT2 = dt_wrt_dL_wrt_dotT2 - dL_wrt_T2;
EQNmotionT2 = expand(EQNmotionT2);
EQNmotionT2 = simplify(EQNmotionT2)