close all;
clear;
clc;

VCase = 0; %gravity =0, gravity and stiffness = 1 TODO: FIX,
NumLinks = 3;

%%%%%%% This version applies pure Lagrangian Mechanics, no attempt to
%%%%%%% define DH Frames.

%% Define World
g = 0;


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

% Damping
b = 10;

% Evaluate the material properties of the points
[m1,I1,k1,~] = MaterialProperties(radius,rho,L1,E);
[m2,I2,k2,~] = MaterialProperties(radius,rho,L2,E);
[m3,I3,k3,~] = MaterialProperties(radius,rho,L3,E);
[~,~,k,~] = MaterialProperties(radius,rho,L1+L2,E);

% Define Magnitude of Magentic moments
mag1 = 1;                             % magnitude of magnetic moment of point 1
mag2 = 1;                             % magnitude of magnetic moment of point 2
mag3 = 1;                             % magnitude of magnetic moment of point 3

%% Define Magnetic Stuff
mu0 = 4*pi*1e-7; % Permeability of free space
B = [25,0,0]; % Magnetic field strength in Tesla (25 mT) in +z direction

% Create the magnetic field
xLow = -0.05;   yLow = -0.05;   zLow = -0.1;
xHigh = 0.05;   yHigh = 0.05;   zHigh = 0.1;
xIncr = 0.005;  yIncr = 0.005;  zIncr = 0.005; %steps size

[x, y, z] = meshgrid(xLow:xIncr:xHigh, yLow:yIncr:yHigh, zLow:zIncr:zHigh);

BxCoil = B(1)*ones(size(x));
ByCoil = B(2)*ones(size(y));
BzCoil = B(3)*ones(size(z));

threshold = 0.001; %Singularity threshold

%% Define initial Starting Conditions

%angles
Theta1 = deg2rad(0);
Theta2 = deg2rad(0);
Theta3 = deg2rad(0);

%angular velocities
DotTheta1 = 0;
DotTheta2 = 0;
DotTheta3 = 0;

% TODO replace these when move to 3D
y0 = 0; y1 = 0; y2 = 0; y3 = 0;

% Obtain Intial Locations TODO Add in y locations
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
    plot3([x0,x1,x2],[y0,y1,y2],[z0,z1,z2],...
        'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1);
    case 3  
    plot3([x0,x1,x2,x3],[y0,y1,y2,y3],[z0,z1,z2,z3],...
        'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1);
    otherwise
end
drawnow;

%% Evaluate Simulation
for t = 1:totalTime


    %% Evaluate Magnetics of each mass

    switch NumLinks

        case 2 % Double Pendulum

            Point1 = [x1,y1,z1];
            Point2 = [x2,y2,z2];

            % Get vector directions of the moments
            mag1Dir = getDirection([x0,y0,z0], Point1);
            mag2Dir = getDirection(Point1, Point2);

            % Calculate Magnetic Moments
            m1_Vector = mag1 * mag1Dir; % mag1 magnitude along Link1's direction
            m2_Vector = mag2 * mag2Dir; % mag2 magnitude along Link2's direction

            % Find meshgrid indices for each dipole
            [idx1, idy1, idz1] = findClosestGridPoint(x, y, z, Point1);
            [idx2, idy2, idz2] = findClosestGridPoint(x, y, z, Point2);

            % Calculate field components for the dipoles
            [Bx1,By1,Bz1] = evaluateField(x,y,z,Point1,mu0,m1_Vector,threshold);
            [Bx2,By2,Bz2] = evaluateField(x,y,z,Point2,mu0,m2_Vector,threshold);

            % Sum the magnetic fields from all sources
            Bx_total = BxCoil + Bx1 + Bx2;
            By_total = ByCoil + By1 + By2;
            Bz_total = BzCoil + Bz1 + Bz2;

            % Calculate magnetic forces
            F1 = f_getForce(Bx_total, By_total, Bz_total, idx1, idy1, idz1, m1_Vector, xIncr, yIncr, zIncr);
            F2 = f_getForce(Bx_total, By_total, Bz_total, idx2, idy2, idz2, m2_Vector, xIncr, yIncr, zIncr);

            % Calculate magnetic torques
            T1 = f_getTorque(BxCoil,ByCoil,BzCoil,idx1, idy1, idz1, m1_Vector);
            T2 = f_getTorque(BxCoil,ByCoil,BzCoil,idx2, idy2, idz2, m2_Vector);
     
        case 3 % Triple Pendulum

            Point1 = [x1,y1,z1];
            Point2 = [x2,y2,z2];
            Point3 = [x3,y3,z3];

            % Get vector directions of the moments
            mag1Dir = getDirection([x0,y0,z0], Point1);
            mag2Dir = getDirection(Point1, Point2);
            mag3Dir = getDirection(Point2, Point3);

            % Calculate Magnetic Moments
            m1_Vector = mag1 * mag1Dir; % mag1 magnitude along Link1's direction
            m2_Vector = mag2 * mag2Dir; % mag2 magnitude along Link2's direction
            m3_Vector = mag3 * mag3Dir; % mag2 magnitude along Link2's direction

            % Find meshgrid indices for each dipole
            [idx1, idy1, idz1] = findClosestGridPoint(x, y, z, Point1);
            [idx2, idy2, idz2] = findClosestGridPoint(x, y, z, Point2);
            [idx3, idy3, idz3] = findClosestGridPoint(x, y, z, Point3);

            % Calculate field components for the dipoles
            [Bx1,By1,Bz1] = evaluateField(x,y,z,Point1,mu0,m1_Vector,threshold);
            [Bx2,By2,Bz2] = evaluateField(x,y,z,Point2,mu0,m2_Vector,threshold);
            [Bx3,By3,Bz3] = evaluateField(x,y,z,Point3,mu0,m3_Vector,threshold);
        
            % Sum the magnetic fields from all sources
            Bx_total = BxCoil + Bx1 + Bx2 + Bx3;
            By_total = ByCoil + By1 + By2 + By3;
            Bz_total = BzCoil + Bz1 + Bz2 + Bz3;

            % Calculate magnetic forces
            F1 = f_getForce(Bx_total, By_total, Bz_total, idx1, idy1, idz1, m1_Vector, xIncr, yIncr, zIncr);
            F2 = f_getForce(Bx_total, By_total, Bz_total, idx2, idy2, idz2, m2_Vector, xIncr, yIncr, zIncr);
            F3 = f_getForce(Bx_total, By_total, Bz_total, idx3, idy3, idz3, m3_Vector, xIncr, yIncr, zIncr);

            % Calculate magnetic torques
            T1 = f_getTorque(BxCoil,ByCoil,BzCoil,idx1, idy1, idz1, m1_Vector);
            T2 = f_getTorque(BxCoil,ByCoil,BzCoil,idx2, idy2, idz2, m2_Vector);
            T3 = f_getTorque(BxCoil,ByCoil,BzCoil,idx3, idy3, idz3, m3_Vector);

        otherwise 

            Point1 = [x1,y1,z1];
            Point2 = [x2,y2,z2];

            % Get vector directions of the moments
            mag1Dir = getDirection([x0,y0,z0], Point1);
            mag2Dir = getDirection(Point1, Point2);

            % Calculate Magnetic Moments
            m1_Vector = mag1 * mag1Dir; % mag1 magnitude along Link1's direction
            m2_Vector = mag2 * mag2Dir; % mag2 magnitude along Link2's direction

            % Find meshgrid indices for each dipole
            [idx1, idy1, idz1] = findClosestGridPoint(x, y, z, Point1);
            [idx2, idy2, idz2] = findClosestGridPoint(x, y, z, Point2);

            % Calculate field components for the dipoles
            [Bx1,By1,Bz1] = evaluateField(x,y,z,Point1,mu0,m1_Vector,threshold);
            [Bx2,By2,Bz2] = evaluateField(x,y,z,Point2,mu0,m2_Vector,threshold);

            % Sum the magnetic fields from all sources
            Bx_total = BxCoil + Bx1 + Bx2;
            By_total = ByCoil + By1 + By2;
            Bz_total = BzCoil + Bz1 + Bz2;

            % Calculate magnetic forces
            F1 = f_getForce(Bx_total, By_total, Bz_total, idx1, idy1, idz1, m1_Vector, xIncr, yIncr, zIncr);
            F2 = f_getForce(Bx_total, By_total, Bz_total, idx2, idy2, idz2, m2_Vector, xIncr, yIncr, zIncr);

            % Calculate magnetic torques
            T1 = f_getTorque(BxCoil,ByCoil,BzCoil,idx1, idy1, idz1, m1_Vector);
            T2 = f_getTorque(BxCoil,ByCoil,BzCoil,idx2, idy2, idz2, m2_Vector);         
    end

    %% Evaluate Dynamics
    
    switch NumLinks
        case 2
            [Theta1,Theta2,~,DotTheta1,DotTheta2,~] = RungeKutte(g,dt,m1,L1,m2,L2,0,0,b,Theta1,Theta2,0,DotTheta1,DotTheta2,0,F1,F2,0,T1,T2,0,NumLinks);
        case 3                                                                
            [Theta1,Theta2,Theta3,DotTheta1,DotTheta2,DotTheta3] = RungeKutte(g,dt,m1,L1,m2,L2,m3,L3,b,Theta1,Theta2,Theta3,DotTheta1,DotTheta2,DotTheta3,F1,F2,F3,T1,T2,T3,NumLinks);
        otherwise
    end
    
    % Determine locations of masses
    [x0,x1,x2,x3,z0,z1,z2,z3] = getPosition(Theta1,Theta2,Theta3,L1,L2,L3);
    
    
    %% Plot the current state
    cla(ax); % Clear the current figure
    
    quiver3(x, y, z, Bx_total, By_total, Bz_total, 0.25, 'b','LineWidth', 0.25);
    hold on
    switch NumLinks
        case 2
        plot3([x0,x1,x2],[y0,y1,y2],[z0,z1,z2],...
            'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1);
        case 3  
        plot3([x0,x1,x2,x3],[y0,y1,y2,y3],[z0,z1,z2,z3],...
            'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1);
        otherwise
    end
        
    drawnow;
    pause(0.1);

end


%% Functions


% Lagrangian gravity only 2 bob
function [DotTheta1, DotTheta2, g1, g2] = LagrangianRHS(m1,L1,m2,L2,g,b,F1,F2,T1,T2,RK)

% This function computes the right hand side of the Euler-Lagrange
% Equations This function uses gravity only for potential energy
% and only works for a double pendulum

    Theta1 = RK(1);
    Theta2 = RK(2);
    DotTheta1 = RK(3);
    DotTheta2 = RK(4);
    
    a1 = (L2/L1)*(m2/(m1+m2))*cos(Theta1-Theta2);
    a2 = (L1/L2)*cos(Theta1-Theta2);
    
    f1 = -(L2/L1)*(m2/(m1+m2))*sin(Theta1-Theta2)*DotTheta2^2 - (g/L1)*sin(Theta1) - b*DotTheta1 + T1(2) + F1(1)*L1*cos(Theta1) + F1(3)*L1*sin(Theta1);
    f2 = (L1/L2)*sin(Theta1-Theta2)*DotTheta1^2 - (g/L2)*sin(Theta2) - b*DotTheta2 + T2(2) + F2(1)*L2*cos(Theta2-Theta1) + F2(3)*L1*sin(Theta2);
    
    g1 = (f1 - (a1*f2))/(1-(a1*a2));
    g2 = (f2 - (a2*f1))/(1-(a1*a2));

end


% Lagrangian gravity and damping only 3 bob                             
function [DotTheta1, DotTheta2, DotTheta3, g1, g2, g3] = LagrangianRHS3(m1,L1,m2,L2,m3,L3,g,b,F1,F2,F3,T1,T2,T3,RK)

% This function computes the right hand side of the Euler-Lagrange
% Equations.


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

    %T1(1) = 0; T2(1) = 0; T3(1) = 0;
    F1 = [0;0;0]; F2 = [0;0;0]; F3 = [0;0;0]; 
    
    f1 = -(g/L1)*sin(Theta1) - T1(2) - F1(1)*L1*cos(Theta1) - F1(3)*L1*sin(Theta1) - (L2/L1)*((m2+m3)/(m1+m2+m3))*sin(Theta1-Theta2)*DotTheta2^2 - (L3/L1)*((m3)/(m1+m2+m3))*sin(Theta1-Theta3)*DotTheta3^2 - b*DotTheta1;
    f2 = -(g/L2)*sin(Theta2) - T2(2) - F2(1)*L2*cos(Theta2) - F2(3)*L1*sin(Theta2) + (L1/L2)*sin(Theta1-Theta2)*DotTheta2^2 - (L3/L2)*((m3)/(m2+m3))*sin(Theta2-Theta3)*DotTheta3^2 - b*DotTheta2;
    f3 = -(g/L3)*sin(Theta3) - T3(2) - F3(1)*L1*cos(Theta3) - F3(3)*L1*sin(Theta3) + (L1/L3)*sin(Theta1-Theta3)*DotTheta1^2 + (L2/L3)*sin(Theta2-Theta3)*DotTheta2^2 - b*DotTheta3;

    F = [f1;f2;f3];

    g = SVDSolve(A, F);

    g1 = g(1);
    g2 = g(2);
    g3 = g(3);

end

                                                                           
function [Theta1,Theta2,Theta3,DotTheta1,DotTheta2,DotTheta3] = RungeKutte(g,dt,m1,L1,m2,L2,m3,L3,b,Theta1,Theta2,Theta3,DotTheta1,DotTheta2,DotTheta3,F1,F2,F3,T1,T2,T3,NumLinks)

% This function applies the Runge-kutte method to obtain the new angles and
% angular velocities. 
    
    switch NumLinks
        case 2 %Double Pendulum
            RK = [Theta1,Theta2,DotTheta1,DotTheta2];
            % Obtain Runge-kutte constants
            [k1(1),k1(2),k1(3),k1(4)] = LagrangianRHS(m1,L1,m2,L2,g,b,F1,F2,T1,T2,RK);
            [k2(1),k2(2),k2(3),k2(4)] = LagrangianRHS(m1,L1,m2,L2,g,b,F1,F2,T1,T2,(RK + dt*k1/2));
            [k3(1),k3(2),k3(3),k3(4)] = LagrangianRHS(m1,L1,m2,L2,g,b,F1,F2,T1,T2,(RK + dt*k2/2));
            [k4(1),k4(2),k4(3),k4(4)] = LagrangianRHS(m1,L1,m2,L2,g,b,F1,F2,T1,T2,(RK + dt*k3));

        case 3 %Triple Pendulum
            RK = [Theta1,Theta2,Theta3,DotTheta1,DotTheta2,DotTheta3];
            % Obtain Runge-kutte constants
            [k1(1),k1(2),k1(3),k1(4),k1(5),k1(6)] = LagrangianRHS3(m1,L1,m2,L2,m3,L3,g,b,F1,F2,F3,T1,T2,T3,RK);
            [k2(1),k2(2),k2(3),k2(4),k2(5),k2(6)] = LagrangianRHS3(m1,L1,m2,L2,m3,L3,g,b,F1,F2,F3,T1,T2,T3,(RK + dt*k1/2));
            [k3(1),k3(2),k3(3),k3(4),k3(5),k3(6)] = LagrangianRHS3(m1,L1,m2,L2,m3,L3,g,b,F1,F2,F3,T1,T2,T3,(RK + dt*k2/2));
            [k4(1),k4(2),k4(3),k4(4),k4(5),k4(6)] = LagrangianRHS3(m1,L1,m2,L2,m3,L3,g,b,F1,F2,F3,T1,T2,T3,(RK + dt*k3));

        otherwise %Double Pendulum
            RK = [Theta1,Theta2,DotTheta1,DotTheta2];
            % Obtain Runge-kutte constants
            [k1(1),k1(2),k1(3),k1(4)] = LagrangianRHS(m1,L1,m2,L2,g,b,F1,F2,T1,T2,RK);
            [k2(1),k2(2),k2(3),k2(4)] = LagrangianRHS(m1,L1,m2,L2,g,b,F1,F2,T1,T2,(RK + dt*k1/2));
            [k3(1),k3(2),k3(3),k3(4)] = LagrangianRHS(m1,L1,m2,L2,g,b,F1,F2,T1,T2,(RK + dt*k2/2));
            [k4(1),k4(2),k4(3),k4(4)] = LagrangianRHS(m1,L1,m2,L2,g,b,F1,F2,T1,T2,(RK + dt*k3));
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
    
    % Evaluate the world space position of each mass.

    % TODO Add in y locations

    x0=0;
    z0=0;

    x1 = L1*sin(Theta1);
    x2 = x1 + L2*sin(Theta2);
    x3 = x2 + L3*sin(Theta3);

    z1 = -L1*cos(Theta1);
    z2 = z1 - L2*cos(Theta2);
    z3 = z2 - L2*cos(Theta3);
end

function g = SVDSolve(A, F)
    % Use singluar value decomposition to solve g = A^-1 * F
    
    % Compute the SVD of A
    [U, S, V] = svd(A);

    % Find the reciprocal of the non-zero elements in S
    tol = max(size(A)) * eps(norm(S, 'fro'));
    S_inv = diag(1 ./ diag(S));
    S_inv(diag(S) < tol) = 0; % Set small values to zero to avoid division by small numbers

    % Compute the pseudo-inverse of A
    A_pinv = V * S_inv * U';

    % Solve for x
    g = A_pinv * F;
end

function dir = getDirection(P1, P2)
    dir = (P2 - P1) / norm(P2 - P1); % Direction from Point1 to Point2;
end

function force = f_getForce(Bx, By, Bz, x_index, y_index, z_index, m, dx, dy, dz)

    % Ensure indices are within valid range
    x_index_minus = max(x_index-1, 1);
    x_index_plus = min(x_index+1, size(Bx,1));
    y_index_minus = max(y_index-1, 1);
    y_index_plus = min(y_index+1, size(By,2));
    z_index_minus = max(z_index-1, 1);
    z_index_plus = min(z_index+1, size(Bz,3));
    
    % Calculate the magnetic field at the dipole's position
    B = [Bx(x_index, y_index, z_index), By(x_index, y_index, z_index), Bz(x_index, y_index, z_index)];

    % Approximate the gradient of the magnetic field using boundary-safe indices
    dBx_dx = (Bx(x_index_plus, y_index, z_index) - Bx(x_index_minus, y_index, z_index)) / (2*dx);
    dBx_dy = (Bx(x_index, y_index_plus, z_index) - Bx(x_index, y_index_minus, z_index)) / (2*dy);
    dBx_dz = (Bx(x_index, y_index, z_index_plus) - Bx(x_index, y_index, z_index_minus)) / (2*dz);
    dBy_dy = (By(x_index, y_index_plus, z_index) - By(x_index, y_index_minus, z_index)) / (2*dy);
    dBz_dz = (Bz(x_index, y_index, z_index_plus) - Bz(x_index, y_index, z_index_minus)) / (2*dz);

    %moment matrix
    M = [m(1) m(2) m(3) 0 0;
         0 m(1) 0 m(2) m(3);
         -m(1) 0 m(1) -m(3) m(2)];

    Grad = [dBx_dx;dBx_dy;dBx_dz;dBy_dy;dBz_dz];

    % Calculate the force
    force = M*Grad;

    % remove singularities
    if isnan(force(1))
        force(1) = 0;
    end
    
    if isnan(force(2))
        force(2) = 0;
    end
    
    if isnan(force(3))
        force(3) = 0;
    end
end

function torque = f_getTorque(Bx, By, Bz, x_index, y_index, z_index, m)
    %This function returns the torques at a dipole wrt another dipole
    
    % Calculate the magnetic field at dipole due to another dipole
    B = [Bx(x_index,y_index,z_index), By(x_index,y_index,z_index), Bz(x_index,y_index,z_index)];
    torque = cross(m, B);
    
    % remove singularities
    if isnan(torque(1))
        torque(1) = 0;
    end
    
    if isnan(torque(2))
        torque(2) = 0;
    end
    
    if isnan(torque(3))
        torque(3) = 0;
    end

end

function [Bx,By,Bz] = evaluateField(x,y,z,Point,mu0,m_vector,threshold)
    %evaluate magnetic field contribution of dipole
    x1 = x - Point(1);
    y1 = y - Point(2);
    z1 = z - Point(3);
    r1 = sqrt(x1.^2 + y1.^2 + z1.^2);
    rx1 = x1./r1; ry1 = y1./r1; rz1 = z1./r1;
    
    Bx = mu0/(4*pi) * (3*(m_vector(1)*rx1 + m_vector(2)*ry1 + m_vector(3)*rz1).*rx1 - m_vector(1))./r1.^3;
    By = mu0/(4*pi) * (3*(m_vector(1)*rx1 + m_vector(2)*ry1 + m_vector(3)*rz1).*ry1 - m_vector(2))./r1.^3;
    Bz = mu0/(4*pi) * (3*(m_vector(1)*rx1 + m_vector(2)*ry1 + m_vector(3)*rz1).*rz1 - m_vector(3))./r1.^3;

    % Remove singularities for all dipoles
    Bx(r1<threshold) = 0; By(r1<threshold) = 0; Bz(r1<threshold) = 0;
end