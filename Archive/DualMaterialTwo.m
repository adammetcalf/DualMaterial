close all;
clear;
clc;

% Pure NdFEB Script

%% Define the world

WorldOrigin = [0,0];
Gravity = [0,-9.81];                  % m/s^2

%% Define the rigid body points

Point1 = [0.01,0];                   % initial condition for point 1
Point2 = [0.02,0];                   % initial condition for point 2
Point3 = [0.03,0];                   % initial condition for point 3
Point4 = [0.04,0];                   % initial condition for point 3
Mass = 0.001;                         % kg


m1 = 0.05;                            % magnitude of magnetic moment of point 1
m2 = 0.05;                            % magnitude of magnetic moment of point 2
m3 = 0.05;                            % magnitude of magnetic moment of point 3
m4 = 0.05;                            % magnitude of magnetic moment of point 3

% Initialise Velocities               % m/s
Velocity1 = [0,0];
Velocity2 = [0,0];
Velocity3 = [0,0];
Velocity4 = [0,0];

%% Define Damping Factor
DampingFactor = 0.01; 

%% Define rigid body links according to the initial conditions.

% Define Links as vectors
Link1 = Point1 - WorldOrigin;
Link2 = Point2 - Point1;
Link3 = Point3 - Point2;
Link4 = Point4 - Point3;

% Assuming Links are rigid, their lengths remain constant
Link1_Length = norm(Link1);
Link2_Length = norm(Link2);
Link3_Length = norm(Link3);
Link4_Length = norm(Link4);

%% Time Setup
dt = 0.01;
totalTime = 2; 
Iterations = totalTime/dt;

%% Define Magnetic Stuff
mu0 = 4*pi*1e-7;                       % Permeability of free space
B = [0,25e-3];                         % Magnetic field strength in Tesla (25 mT) in +y direction

% Create the magnetic field
xLow = -0.05;   yLow = -0.05;   
xHigh = 0.05;   yHigh = 0.05;   
xIncr = 0.0002;  yIncr = 0.0002; %steps size

[x, y] = meshgrid(xLow:xIncr:xHigh, yLow:yIncr:yHigh);

BxCoil = B(1)*ones(size(x));
ByCoil = B(2)*ones(size(y));


threshold = 0.001; %Singularity threshold

% Create a figure for the animation
hFig = figure;

% Set the figure to fullscreen
set(hFig, 'units', 'normalized', 'outerposition', [0 0 1 1]);

% Set up the axes for the 3D plot
ax = axes('Parent', hFig);
hold(ax, 'on');
axis(ax, 'equal');
axis(ax,[-0.05 0.05 -0.05 0.05]);  
grid(ax, 'on');
xlabel(ax, 'X');
ylabel(ax, 'Y');
title(ax, 'Dual Material Sim');


%% Simulation Loop
for step = 1:Iterations

    m1_dir = getDirection(WorldOrigin, Point1);
    m2_dir = getDirection(Point1, Point2);
    m3_dir = getDirection(Point2, Point3);
    m4_dir = getDirection(Point3, Point4);


    % Calculate Magnetic Moments
    m1_Vector = m1 * m1_dir; % m1 magnitude along Link1's direction
    m2_Vector = m2 * m2_dir; % m2 magnitude along Link1's direction
    m3_Vector = m3 * m3_dir; % m3 magnitude along Link1's direction
    m4_Vector = m4 * m4_dir; % m4 magnitude along Link1's direction

    % Calculate field components for the dipoles
    [Bx1,By1] = evaluateField(x,y,Point1,mu0,m1_Vector,threshold);
    [Bx2,By2] = evaluateField(x,y,Point2,mu0,m2_Vector,threshold);
    [Bx3,By3] = evaluateField(x,y,Point3,mu0,m3_Vector,threshold);
    [Bx4,By4] = evaluateField(x,y,Point4,mu0,m4_Vector,threshold);

    % Sum the magnetic fields from both dipoles
    Bx_total = BxCoil + Bx1 + Bx2 + Bx3 + Bx4;
    By_total = ByCoil + By1 + By2 + By3 + By4;

    % Find meshgrid indices for each dipole
    [idx1, idy1] = findClosestGridPoint(x, y, Point1);
    [idx2, idy2] = findClosestGridPoint(x, y, Point2);
    [idx3, idy3] = findClosestGridPoint(x, y, Point3);
    [idx4, idy4] = findClosestGridPoint(x, y, Point4);

    % Calculate magnetic forces
    F1 = f_getForce(Bx_total, By_total, idx1, idy1, m1_Vector, xIncr, yIncr);
    F2 = f_getForce(Bx_total, By_total, idx2, idy2, m2_Vector, xIncr, yIncr);
    F3 = f_getForce(Bx_total, By_total, idx3, idy3, m3_Vector, xIncr, yIncr);
    F4 = f_getForce(Bx_total, By_total, idx4, idy4, m4_Vector, xIncr, yIncr);


    % Calculate accelerations due to magnetic forces
    acceleration1 = F1 / Mass + Gravity;
    acceleration2 = F2 / Mass + Gravity;
    acceleration3 = F3 / Mass + Gravity;
    acceleration4 = F4 / Mass + Gravity;

    Velocity1 = Velocity1 + (acceleration1 - DampingFactor * Velocity1) * dt;
    Velocity2 = Velocity2 + (acceleration1 - DampingFactor * Velocity2) * dt;
    Velocity3 = Velocity3 + (acceleration1 - DampingFactor * Velocity3) * dt;
    Velocity4 = Velocity4 + (acceleration1 - DampingFactor * Velocity4) * dt;

    %% Update positions based on velocities
    Point1 = Point1 + Velocity1 * dt;
    Point2 = Point2 + Velocity2 * dt;
    Point3 = Point3 + Velocity3 * dt;
    Point4 = Point4 + Velocity4 * dt;

    % Enforce link constraints and get magentic moment direction
    Point1 = enforceLinkConstraint(WorldOrigin, Point1, Link1_Length);
    Point2 = enforceLinkConstraint(Point1, Point2, Link2_Length);
    Point3 = enforceLinkConstraint(Point2, Point3, Link3_Length);
    Point4 = enforceLinkConstraint(Point3, Point4, Link4_Length);

    %% Plot the current state
    cla(ax); % Clear the current figure

    plot([WorldOrigin(1), Point1(1), Point2(1), Point3(1), Point4(1)],...
          [WorldOrigin(2), Point1(2), Point2(2), Point3(2), Point4(2)],...
          'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1);
    hold on
    %plot Magnetic field
    quiver(x, y, Bx_total, By_total, 2, 'b','LineWidth', 0.5);
    hold on;
    
    drawnow;

    % Capture frames for a video
    frames(step) = getframe(hFig);

    pause(dt);
end

% Save the animation as a video file
video = VideoWriter('DualMaterial2.avi');
video.FrameRate = 10;
open(video);
writeVideo(video, frames);
close(video);


%% Functions

function dir = getDirection(P1, P2)
    dir = (P2 - P1) / norm(P2 - P1); % Direction from Point1 to Point2;
end


function [Bx,By] = evaluateField(x,y,Point,mu0,m_vector,threshold)
    %evaluate magnetic field contribution of dipole
    x1 = x - Point(1);
    y1 = y - Point(2);
    r1 = sqrt(x1.^2 + y1.^2);
    rx1 = x1./r1; ry1 = y1./r1;
    
    Bx = mu0/(4*pi) * (3*(m_vector(1)*rx1 + m_vector(2)*ry1).*rx1 - m_vector(1))./r1.^3;
    By = mu0/(4*pi) * (3*(m_vector(1)*rx1 + m_vector(2)*ry1).*ry1 - m_vector(2))./r1.^3;
    
    % Remove singularities for all dipoles
    Bx(r1<threshold) = 0; By(r1<threshold) = 0;
end

% Function to find the closest grid point
function [idx, idy] = findClosestGridPoint(x, y, p)
    [~, idx] = min(abs(x(1,:,1) - p(1)));
    [~, idy] = min(abs(y(:,1,1) - p(2)));
end

function force = f_getForce(Bx, By, x_index, y_index, m, dx, dy)
    % Ensure indices are within valid range
    x_index_minus = max(x_index-1, 1);
    x_index_plus = min(x_index+1, size(Bx,1));
    y_index_minus = max(y_index-1, 1);
    y_index_plus = min(y_index+1, size(By,2));
    
    % Calculate the magnetic field at the dipole's position
    B = [Bx(x_index, y_index), By(x_index, y_index)];

    % Approximate the gradient of the magnetic field using boundary-safe indices
    dBx_dx = (Bx(x_index_plus, y_index) - Bx(x_index_minus, y_index)) / (2*dx);
    dBy_dy = (By(x_index, y_index_plus) - By(x_index, y_index_minus)) / (2*dy);


    % Calculate the force
    force = [m(1) * dBx_dx, m(2) * dBy_dy];

    % remove singularities
    if isnan(force(1))
        force(1) = 0;
    end
    
    if isnan(force(2))
        force(2) = 0;
    end
    
end

function Pos2 = enforceLinkConstraint(P1, P2, L2)
    % Enforce the distance between Point1 and Point2
    dir = (P2 - P1) / norm(P2 - P1); % Direction from Point1 to Point2
    Pos2 = P1 + dir * L2; % Adjust Point2 position correctly
end