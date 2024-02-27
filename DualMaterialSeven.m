close all;
clear;
clc;

%%%%%%% This version applies pure Lagrangian Mechanics, no attempt to
%%%%%%% define DH Frames.

%% Define World
g = 9.81;

%% Define material properties

% Lengths
L1 = 0.01; %Length 1 in meters
L2 = 0.01; %Length 2 in meters

% Link Radii (m)
radius = 1e-03;

% EcoFlex 0030 Density (kg m^-3)
rho = 1070; 

% Youngs Modulus Ecoflex 0030 (125 kPa)
E = 125000;

% Evaluate the material properties of the points
[m1,I1,k1,~] = MaterialProperties(radius,rho,L1,E);
[m2,I2,k2,~] = MaterialProperties(radius,rho,L2,E);

%% Define initial Starting Conditions

%angles
Theta1 = deg2rad(90);
Theta2 = deg2rad(90);

%angular velocities
DotTheta1 = 0;
DotTheta2 = 0;

[x0,x1,x2,z0,z1,z2] = getPosition(Theta1,Theta2,L1,L2);

%% Check Variables so far
% Runge-kutta method requires m1, L1, Theta1, thetaDot1, m2, L2, Theta2, thetaDot2 

%% Time Setup
dt = 0.01;
totalTime = 100; 
Iterations = totalTime/dt;

%% Animation Setup

% Create a figure for the animation
hFig = figure;

% Set the figure to fullscreen
set(hFig, 'units', 'normalized', 'outerposition', [0 0 1 1]);

% Set up the axes for the 3D plot
ax = axes('Parent', hFig);
hold(ax, 'on');
axis(ax, 'equal');
axis(ax,[-0.05 0.05 -0.05 0.05 -0.05 0.05]);  
grid(ax, 'on');
xlabel(ax, 'X');
ylabel(ax, 'Y');
zlabel(ax, 'Z');
title(ax, 'Dual Material Sim');
view(22,7);

cla(ax); % Clear the current figure

plot3([x0,x1,x2],[0,0,0],[z0,z1,z2],...
          'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1);
hold on   
drawnow;

for t = 1:totalTime

    %% Plot the current state
   
    [Theta1,Theta2,DotTheta1,DotTheta2] = RungeKutte(g,dt,m1,L1,m2,L2,Theta1,Theta2,DotTheta1,DotTheta2);
    [x0,x1,x2,z0,z1,z2] = getPosition(Theta1,Theta2,L1,L2);

    cla(ax); % Clear the current figure

    plot3([x0,x1,x2],[0,0,0],[z0,z1,z2],...
          'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1);
    hold on   
    drawnow;
    pause(0.1);

end






function [DotTheta1, DotTheta2, g1, g2] = LagrangianRHS(m1,L1,m2,L2,g,RK)

% This function computes the right hand side of the Euler-Lagrange
% Equations

    Theta1 = RK(1);
    Theta2 = RK(2);
    DotTheta1 = RK(3);
    DotTheta2 = RK(4);
    
    a1 = (L2/L1)*(m2/(m1+m2))*cos(Theta1-Theta2);
    a2 = (L1/L2)*cos(Theta1-Theta2);
    
    f1 = -(L2/L1)*(m2/(m1+m2))*sin(Theta1-Theta2)*DotTheta2^2 - (g/L1)*sin(Theta1);
    f2 = (L1/L2)*sin(Theta1-Theta2)*DotTheta1^2 - (g/L2)*sin(Theta2);
    
    g1 = (f1 - (a1*f2))/(1-(a1*a2));
    g2 = (f2 - (a2*f1))/(1-(a1*a2));

end

function [Theta1,Theta2,DotTheta1,DotTheta2] = RungeKutte(g,dt,m1,L1,m2,L2,Theta1,Theta2,DotTheta1,DotTheta2)

% This function applies the Runge-kutte method to obtain the new angles and
% angular velocities

    RK = [Theta1,Theta2,DotTheta1,DotTheta2];
    
    % Obtain Runge-kutte constants
    [k1(1),k1(2),k1(3),k1(4)] = LagrangianRHS(m1,L1,m2,L2,g,RK);
    [k2(1),k2(2),k2(3),k2(4)] = LagrangianRHS(m1,L1,m2,L2,g,(RK + dt*k1/2));
    [k3(1),k3(2),k3(3),k3(4)] = LagrangianRHS(m1,L1,m2,L2,g,(RK + dt*k2/2));
    [k4(1),k4(2),k4(3),k4(4)] = LagrangianRHS(m1,L1,m2,L2,g,(RK + dt*k3));
    
    % Compute the RK4 right hand side
    R = (1/6)*dt*(k1 + 2*k2 + 2*k3 + k4);
    
    Theta1 = Theta1 + R(1);
    Theta2 = Theta2 + R(2);
    DotTheta1 = DotTheta1 + R(3);
    DotTheta2 = DotTheta2 + R(4);

end

function [x0,x1,x2,z0,z1,z2] = getPosition(Theta1,Theta2,L1,L2)

    x0=0;
    z0=0;

    x1 = L1*sin(Theta1);
    x2 = x1 + L2*sin(Theta2);

    z1 = -L1*cos(Theta1);
    z2 = z1 - L2*cos(Theta2);

end