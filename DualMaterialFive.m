close all;
clear;
clc;

% Pure NdFEB Script

%% General DH Frame Transformation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%            Not used, just for documentation            %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Frame = [cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha) a*cos(theta)
%           sin(theta) cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta)
%           0 sin(alpha) cos(alpha) d
%           0 0 0 1];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Jacobean Matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%            Not used, just for documentation            %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Two uses for the Jacobean matrix. 
% 1. Tau = J^T F
% 2. X_dot = J Theta_dot


%% Define Rotation 
% These may be used to adjust the starting position of the tentacle.

% Note Bene, NO LONGER USED, Just here for documentation.

Beta = deg2rad(90);

rotX = [1 0 0 0
        0 cos(Beta) -sin(Beta) 0
        0 sin(Beta) cos(Beta) 0
        0 0 0 1];

rotY = [cos(Beta) 0 sin(Beta) 0
        0 1 0 0
        -sin(Beta) 0 cos(Beta) 0
        0 0 0 1];

rotZ = [cos(Beta) -sin(Beta) 0 0
        sin(Beta) cos(Beta) 0 0
        0 0 1 0
        0 0 0 1];

NoRot = [1 0 0 0
         0 1 0 0
         0 0 1 0
         0 0 0 1];

%% Define DH Frames

theta0 = deg2rad(90); %% Initial World frame Z axis orientation defined by this value
alpha0 = deg2rad(180); %% Initial starting position defined by this value (x axis).
alpha1 = deg2rad(0);
alpha2 = deg2rad(00);
alpha3 = deg2rad(0);
alpha4 = deg2rad(0);

l1 = 0.01;
l2 = 0.01;
l3 = 0.01;
l4 = 0.01;



%% Define World

% Define the effects of gravity (m/s^2)
Gravity = [0,0,-9.81];                

% Define startiung position of tentacle by apply transformation.
StartTransform = NoRot; %%Note, startn position is now defined using theta0


%% Define Tentacle physical Characteristics

% Link Radii (m)
radius = 1e-03;

% EcoFlex 0030 Density (kg m^-3)
rho = 1070; 

% Youngs Modulus Ecoflex 0030 (125 kPa)
E = 125000;

% Evaluate the material properties of the points
[Mass1,I1,k1,Length1] = MaterialProperties(radius,rho,l1,E);
[Mass2,I2,k2,Length2] = MaterialProperties(radius,rho,l2,E);
[Mass3,I3,k3,Length3] = MaterialProperties(radius,rho,l3,E);
[Mass4,I4,k4,Length4] = MaterialProperties(radius,rho,l4,E);


%% Define initial conditions

% Translational velocities (m s^-1)
v1 = [0,0,0]; 
v2 = [0,0,0];
v3 = [0,0,0];
v4 = [0,0,0];


%% Magnetic Properties of the tentacle
% TODO Implement magnetic properties
% Define Magnetic Moments (Magnitudes)
m1 = 0.05;                             % magnitude of magnetic moment of point 1
m2 = 0.05;                             % magnitude of magnetic moment of point 2
m3 = 0.05;                             % magnitude of magnetic moment of point 3
m4 = 0.05;                             % magnitude of magnetic moment of point 4

%% Place Tentacle 

% Evaluate Forward Kinematics and Jacobian
[Origin, Point1, Point2, Point3, Point4,J] = Tent5FK(theta0,alpha0,alpha1,alpha2,alpha3,alpha4,StartTransform,l1,l2,l3,l4);


% Extract Position in World Frame
WorldOrigin = Origin(1:3,4);
Point1Pos = Point1(1:3,4);
Point2Pos = Point2(1:3,4);
Point3Pos = Point3(1:3,4);
Point4Pos = Point4(1:3,4);

dt = 0.01; % Time step in seconds
time = 0:dt:2; % Total simulation time from 0 to 2 seconds

for t = 1

    % Calculate the force of gravity on each segment
    F_gravity1 = Gravity.*Mass1;
    F_gravity2 = Gravity.*Mass2;
    F_gravity3 = Gravity.*Mass3;
    F_gravity4 = Gravity.*Mass4;
    
    % Update velocities of each point
    v1 = v1 + (F_gravity1 / Mass1) * dt;
    v2 = v2 + (F_gravity2 / Mass2) * dt;
    v3 = v3 + (F_gravity3 / Mass3) * dt;
    v4 = v4 + (F_gravity4 / Mass4) * dt;
   
    % Evaluate Forward Kinematics and Jacobian
    [Origin, Point1, Point2, Point3, Point4,J] = Tent5FK(theta0,alpha0,alpha1,alpha2,alpha3,alpha4,StartTransform,l1,l2,l3,l4);

    %% Plot Tentacle
    figure(1)
    clf; % Clear the current figure
    subplot(1,2,1)
    plot3([WorldOrigin(1),Point1Pos(1),Point2Pos(1),Point3Pos(1),Point4Pos(1)],...
          [WorldOrigin(2),Point1Pos(2),Point2Pos(2),Point3Pos(2),Point4Pos(2)],...
          [WorldOrigin(3),Point1Pos(3),Point2Pos(3),Point3Pos(3),Point4Pos(3)],...
           'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1)
    hold on
    axis equal;
    axis([-0.05 0.05 -0.05 0.05 -0.05 0.05]);
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    subplot(1,2,2)
    scale = 0.005;
    plot3([WorldOrigin(1),Point1Pos(1),Point2Pos(1),Point3Pos(1),Point4Pos(1)],...
          [WorldOrigin(2),Point1Pos(2),Point2Pos(2),Point3Pos(2),Point4Pos(2)],...
          [WorldOrigin(3),Point1Pos(3),Point2Pos(3),Point3Pos(3),Point4Pos(3)],...
           'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1)
    hold on
    addOrientationArrows(Origin,scale);
    addOrientationArrows(Point1,scale);
    addOrientationArrows(Point2,scale);
    addOrientationArrows(Point3,scale);
    addOrientationArrows(Point4,scale);
    axis equal;
    axis([-0.05 0.05 -0.05 0.05 -0.05 0.05]);
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    % Set the figure to full screen using WindowState
    set(gcf, 'WindowState', 'maximized');
    drawnow;
end

