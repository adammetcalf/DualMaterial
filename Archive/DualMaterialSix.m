close all
clear
clc;

syms theta_1(t) theta_2(t) theta_3(t) L_1 L_2 L_3 m_1 m_2 m_3 g

% Define displacements in cartesian coordinates
x_1 = L_1*sin(theta_1);
y_1 = -L_1*cos(theta_1);
x_2 = x_1 + L_2*sin(theta_2);
y_2 = y_1 - L_2*cos(theta_2);
x_3 = x_2 + L_3*sin(theta_3);
y_3 = y_2 - L_3*cos(theta_3);

% Differentiate displacements to find velocities
vx_1 = diff(x_1);
vy_1 = diff(y_1);
vx_2 = diff(x_2);
vy_2 = diff(y_2);
vx_3 = diff(x_3, t);
vy_3 = diff(y_3, t);

% Differentiate velocities to find accelerations
ax_1 = diff(vx_1);
ay_1 = diff(vy_1);
ax_2 = diff(vx_2);
ay_2 = diff(vy_2);
ax_3 = diff(vx_3, t);
ay_3 = diff(vy_3, t);

% Define rod tensions
syms T_1 T_2 T_3

% Define forces on mass 1
eqx_1 = m_1*ax_1(t) == -T_1*sin(theta_1(t)) + T_2*sin(theta_2(t));
eqy_1 = m_1*ay_1(t) == T_1*cos(theta_1(t)) - T_2*cos(theta_2(t)) - m_1*g;

% Define forces on mass 2
eqx_2 = m_2*ax_2(t) == -T_2*sin(theta_2(t));
eqy_2 = m_2*ay_2(t) == T_2*cos(theta_2(t)) - m_2*g;

% Define forces on mass 3
eqx_3 = m_3*ax_3(t) == -T_3*sin(theta_3(t));
eqy_3 = m_3*ay_3(t) == T_3*cos(theta_3(t)) - m_3*g;

% Evaluate Tension forces
Tension = solve([eqx_1 eqy_1],[T_1 T_2 T_3]);

% Substitute Tension forces into mass 3
eqRed_1 = subs(eqx_3,[T_1 T_2 T_3],[Tension.T_1 Tension.T_2 Tension.T_3]);
eqRed_2 = subs(eqy_3,[T_1 T_2 T_3],[Tension.T_1 Tension.T_2 Tension.T_3]);

% Solve system equations
L_1 = 1;
L_2 = 1;
L_3 = 1;
m_1 = 1;
m_2 = 1;
m_3 = 1;
g = 9.8;

eqn_1 = subs(eqRed_1);
eqn_2 = subs(eqRed_2);

% Convert to 1st order equations (currently they are non linear second
% order)
[V,S] = odeToVectorField(eqn_1,eqn_2);

% Further convert to a MATLAB Function
M = matlabFunction(V,'vars',{'t','Y'});

% Define initial conditions
initCond = [pi/4 0 pi/6 0 pi/4 0];

% Solve for state variables
sols = ode45(M,[0 10],initCond);


% Animate
x_1 = @(t) L_1*sin(deval(sols,t,3));
y_1 = @(t) -L_1*cos(deval(sols,t,3));
x_2 = @(t) L_1*sin(deval(sols,t,3))+L_2*sin(deval(sols,t,1));
y_2 = @(t) -L_1*cos(deval(sols,t,3))-L_2*cos(deval(sols,t,1));
x_3 = @(t) L_1*sin(deval(sols,t,3)) + L_2*sin(deval(sols,t,1)) + L_3*sin(deval(sols,t,5)); % Update index 5 if needed
y_3 = @(t) -L_1*cos(deval(sols,t,3)) - L_2*cos(deval(sols,t,1)) - L_3*cos(deval(sols,t,5)); % Update index 5 if needed

fanimator(@(t) plot(x_1(t),y_1(t),'ro','MarkerSize',m_1*10,'MarkerFaceColor','r'));
axis equal;

hold on;
fanimator(@(t) plot([0 x_1(t)],[0 y_1(t)],'r-'));
fanimator(@(t) plot(x_2(t),y_2(t),'go','MarkerSize',m_2*10,'MarkerFaceColor','g'));
fanimator(@(t) plot([x_1(t) x_2(t)],[y_1(t) y_2(t)],'g-'));

fanimator(@(t) text(-0.3,0.3,"Timer: "+num2str(t,2)));
hold off;

playAnimation

