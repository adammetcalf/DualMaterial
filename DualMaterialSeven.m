close all;
clear;
clc;

VCase = 0; %gravity =0, gravity and stiffness = 1 TODO: FIX,
NumLinks = 3;

%%%%%%% This version applies pure Lagrangian Mechanics, no attempt to
%%%%%%% define DH Frames.

%% Define World
g = 9.81;
b = 100;

%% Define material properties

% Lengths
L1 = 0.01; %Length 1 in meters
L2 = 0.01; %Length 2 in meters
L3 = 0.01; %Length 3 in meters

% Link Radii (m)
radius = 1e-03;

% EcoFlex 0030 Density (kg m^-3)
rho = 1070; 

% Youngs Modulus Ecoflex 0030 (125 kPa)
E = 125000;

% Evaluate the material properties of the points
[m1,I1,k1,~] = MaterialProperties(radius,rho,L1,E);
[m2,I2,k2,~] = MaterialProperties(radius,rho,L2,E);
[m3,I3,k3,~] = MaterialProperties(radius,rho,L3,E);
[~,~,k,~] = MaterialProperties(radius,rho,L1+L2,E);

%% Define initial Starting Conditions

%angles
Theta1 = deg2rad(45);
Theta2 = deg2rad(45);
Theta3 = deg2rad(90);

%angular velocities
DotTheta1 = 0;
DotTheta2 = 0;
DotTheta3 = 0;

[x0,x1,x2,x3,z0,z1,z2,z3] = getPosition(Theta1,Theta2,Theta3,L1,L2,L3);

%% Time Setup
dt = 0.01;
totalTime = 1000; 
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

switch NumLinks
    case 2
    plot3([x0,x1,x2],[0,0,0],[z0,z1,z2],...
        'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1);
    case 3  
    plot3([x0,x1,x2,x3],[0,0,0,0],[z0,z1,z2,z3],...
        'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1);
    otherwise
end
hold on   
drawnow;

for t = 1:totalTime

    %% Plot the current state
   
    switch NumLinks
        case 2
            [Theta1,Theta2,~,DotTheta1,DotTheta2,~] = RungeKutte(g,dt,m1,L1,m2,L2,0,0,Theta1,Theta2,0,DotTheta1,DotTheta2,0,NumLinks);
        case 3  
            [Theta1,Theta2,Theta3,DotTheta1,DotTheta2,DotTheta3] = RungeKutte(g,dt,m1,L1,m2,L2,m3,L3,b,Theta1,Theta2,Theta3,DotTheta1,DotTheta2,DotTheta3,NumLinks);
        otherwise
    end
    

    [x0,x1,x2,x3,z0,z1,z2,z3] = getPosition(Theta1,Theta2,Theta3,L1,L2,L3);

    cla(ax); % Clear the current figure

    switch NumLinks
        case 2
        plot3([x0,x1,x2],[0,0,0],[z0,z1,z2],...
            'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1);
        case 3  
        plot3([x0,x1,x2,x3],[0,0,0,0],[z0,z1,z2,z3],...
            'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1);
        otherwise
    end
    hold on 
    drawnow;
    pause(0.1);

end


%% Functions


% Lagrangian gravity only 2 bob
function [DotTheta1, DotTheta2, g1, g2] = LagrangianRHS(m1,L1,m2,L2,g,RK)

% This function computes the right hand side of the Euler-Lagrange
% Equations This function uses gravity only for potential energy
% and only works for a double pendulum

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


% Lagrangian gravity and Stiffness 2 bob
function [DotTheta1, DotTheta2, g1, g2] = LagrangianRHS2(m1,L1,m2,L2,g,k,RK)

% This function computes the right hand side of the Euler-Lagrange
% Equations This function uses gravity and stiffness for potential energy
% and only works for a double pendulum

    Theta1 = RK(1);
    Theta2 = RK(2);
    DotTheta1 = RK(3);
    DotTheta2 = RK(4);

    % NOTE - Adding stiffness like this is making huge bloody numbers that
    % break the simulation. This is due to the low stiffness of the
    % material being divided by very tiny numbers of mass and length. Note
    % that the Equations were formed using Vstiff =
    % (1/2)*k*(Theta2-Theta1)^2 TODO FIX THIS.
        
    
    a1 = (L2/L1)*(m2/(m1+m2))*cos(Theta1-Theta2);
    a2 = (L1/L2)*cos(Theta1-Theta2);
    
    f1 = -(L2/L1)*(m2/(m1+m2))*sin(Theta1-Theta2)*DotTheta2^2 - (g/L1)*sin(Theta1) - (k/((m1+m2)*L1^2))*(Theta1-Theta2);
    f2 = (L1/L2)*sin(Theta1-Theta2)*DotTheta1^2 - (g/L2)*sin(Theta2) - (k/(m2*L2^2))*(Theta2-Theta1);
    
    g1 = (f1 - (a1*f2))/(1-(a1*a2));
    g2 = (f2 - (a2*f1))/(1-(a1*a2));

end

% Lagrangian gravity only 3 bob
function [DotTheta1, DotTheta2, DotTheta3, g1, g2, g3] = LagrangianRHS3(m1,L1,m2,L2,m3,L3,g,b,RK)

% This function computes the right hand side of the Euler-Lagrange
% Equations This function uses gravity only for potential energy
% and only works for a double pendulum

    Theta1 = RK(1);
    Theta2 = RK(2);
    Theta3 = RK(3);
    DotTheta1 = RK(4);
    DotTheta2 = RK(5);
    DotTheta3 = RK(6); 

    a1 = (L2/L1)*((m2/(m1+m2+m3))+(m2/(m1+m2+m3)))*cos(Theta1-Theta2);
    a2 = (L3/L1)*((m3)/(m1+m2+m3))*cos(Theta1-Theta3);
    a3 = (L1/L2)*cos(Theta1-Theta2);
    a4 = (L3/L2)*((m3)/(m2+m3))*cos(Theta2-Theta3);
    a5 = (L1/L3)*cos(Theta1-Theta3);
    a6 = (L2/L3)*cos(Theta2-Theta3);

    A = [1 a1 a2; a3 1 a4; a5 a6 1];
    %invA = inv(A);

    f1 = -(g/L1)*sin(Theta1) - (L2/L1)*((m2+m3)/(m1+m2+m3))*sin(Theta1-Theta2)*DotTheta2^2 - (L3/L1)*((m3)/(m1+m2+m3))*sin(Theta1-Theta3)*DotTheta3^2 - b*DotTheta1;
    f2 = -(g/L2)*sin(Theta2) + (L1/L2)*sin(Theta1-Theta2)*DotTheta2^2 - (L3/L2)*((m3)/(m2+m3))*sin(Theta2-Theta3)*DotTheta3^2 - b*DotTheta1;
    f3 = -(g/L3)*sin(Theta3) + (L1/L3)*sin(Theta1-Theta3)*DotTheta1^2 + (L2/L3)*sin(Theta2-Theta3)*DotTheta2^2 - b*DotTheta1;

    F = [f1;f2;f3];

    g = A\F;

    g1 = g(1);
    g2 = g(2);
    g3 = g(3);

end


function [Theta1,Theta2,Theta3,DotTheta1,DotTheta2,DotTheta3] = RungeKutte(g,dt,m1,L1,m2,L2,m3,L3,b,Theta1,Theta2,Theta3,DotTheta1,DotTheta2,DotTheta3,NumLinks)

% This function applies the Runge-kutte method to obtain the new angles and
% angular velocities. 

    
    
    switch NumLinks
        case 2 %Double Pendulum
            RK = [Theta1,Theta2,DotTheta1,DotTheta2];
            % Obtain Runge-kutte constants
            [k1(1),k1(2),k1(3),k1(4)] = LagrangianRHS(m1,L1,m2,L2,g,RK);
            [k2(1),k2(2),k2(3),k2(4)] = LagrangianRHS(m1,L1,m2,L2,g,(RK + dt*k1/2));
            [k3(1),k3(2),k3(3),k3(4)] = LagrangianRHS(m1,L1,m2,L2,g,(RK + dt*k2/2));
            [k4(1),k4(2),k4(3),k4(4)] = LagrangianRHS(m1,L1,m2,L2,g,(RK + dt*k3));

        case 3 %Triple Pendulum
            RK = [Theta1,Theta2,Theta3,DotTheta1,DotTheta2,DotTheta3];
            % Obtain Runge-kutte constants
            [k1(1),k1(2),k1(3),k1(4),k1(5),k1(6)] = LagrangianRHS3(m1,L1,m2,L2,m3,L3,g,b,RK);
            [k2(1),k2(2),k2(3),k2(4),k2(5),k2(6)] = LagrangianRHS3(m1,L1,m2,L2,m3,L3,g,b,(RK + dt*k1/2));
            [k3(1),k3(2),k3(3),k3(4),k3(5),k3(6)] = LagrangianRHS3(m1,L1,m2,L2,m3,L3,g,b,(RK + dt*k2/2));
            [k4(1),k4(2),k4(3),k4(4),k4(5),k4(6)] = LagrangianRHS3(m1,L1,m2,L2,m3,L3,g,b,(RK + dt*k3));

        otherwise %Double Pendulum
            RK = [Theta1,Theta2,DotTheta1,DotTheta2];
            % Obtain Runge-kutte constants
            [k1(1),k1(2),k1(3),k1(4)] = LagrangianRHS(m1,L1,m2,L2,g,RK);
            [k2(1),k2(2),k2(3),k2(4)] = LagrangianRHS(m1,L1,m2,L2,g,(RK + dt*k1/2));
            [k3(1),k3(2),k3(3),k3(4)] = LagrangianRHS(m1,L1,m2,L2,g,(RK + dt*k2/2));
            [k4(1),k4(2),k4(3),k4(4)] = LagrangianRHS(m1,L1,m2,L2,g,(RK + dt*k3));
    end
    % Compute the RK4 right hand side
    R = (1/6)*dt*(k1 + 2*k2 + 2*k3 + k4);
    
     switch NumLinks
        case 2 %Double Pendulum
            Theta1 = Theta1 + R(1);
            Theta2 = Theta2 + R(2);
            DotTheta1 = DotTheta1 + R(3);
            DotTheta2 = DotTheta2 + R(4);
    
        case 3 %Triple Pendulum
            Theta1 = Theta1 + R(1);
            Theta2 = Theta2 + R(2);
            Theta3 = Theta3 + R(3);
            DotTheta1 = DotTheta1 + R(4);
            DotTheta2 = DotTheta2 + R(5);
            DotTheta3 = DotTheta3 + R(6);
    
        otherwise %Double Pendulum
            Theta1 = Theta1 + R(1);
            Theta2 = Theta2 + R(2);
            DotTheta1 = DotTheta1 + R(3);
            DotTheta2 = DotTheta2 + R(4);
     end
end

function [x0,x1,x2,x3,z0,z1,z2,z3] = getPosition(Theta1,Theta2,Theta3,L1,L2,L3)

    x0=0;
    z0=0;

    x1 = L1*sin(Theta1);
    x2 = x1 + L2*sin(Theta2);
    x3 = x2 + L3*sin(Theta3);

    z1 = -L1*cos(Theta1);
    z2 = z1 - L2*cos(Theta2);
    z3 = z2 - L2*cos(Theta3);
end