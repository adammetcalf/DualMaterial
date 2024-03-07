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

%% Define Rotation 
% These may be used to adjust the starting position of the tentacle.

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

alpha0 = deg2rad(0); %% Initial starting position defined by this value (x axis).
beta0 = deg2rad(0); %% Initial starting position defined by this value (y axis).

[OriginDH_x, OriginRot_z, OriginDH_y, OriginRot_Negz,...
    Point1DH, Point2DH, Point3DH, Point4DH] = Tent4FK(alpha0,beta0);


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
[Mass1,I1,k1,Length1] = MaterialProperties(radius,rho,Point1DH,E);
[Mass2,I2,k2,Length2] = MaterialProperties(radius,rho,Point2DH,E);
[Mass3,I3,k3,Length3] = MaterialProperties(radius,rho,Point3DH,E);
[Mass4,I4,k4,Length4] = MaterialProperties(radius,rho,Point4DH,E);


% TODO Jacobean?

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

% Evaluate Forward Kinematics
Originx = StartTransform*OriginDH_x;
Originy = Originx*OriginRot_z*OriginDH_y;
Point1 = Originy*OriginRot_Negz*Point1DH;
Point2 = Point1*Point2DH;
Point3 = Point2*Point3DH;
Point4 = Point3*Point4DH;

% Extract Position in World Frame
WorldOrigin = Originy(1:3,4);
Point1Pos = Point1(1:3,4);
Point2Pos = Point2(1:3,4);
Point3Pos = Point3(1:3,4);
Point4Pos = Point4(1:3,4);

dt = 0.01; % Time step in seconds
time = 0:dt:2; % Total simulation time from 0 to 2 seconds

for t = 1:90

    alpha0 = deg2rad(t);
    beta0 = deg2rad(t);

    [OriginDH_x, OriginRot_z, OriginDH_y, OriginRot_Negz,...
    Point1DH, Point2DH, Point3DH, Point4DH] = Tent4FK(alpha0,beta0);

    % Evaluate Forward Kinematics
    Originx = StartTransform*OriginDH_x;
    Originy = Originx*OriginRot_z*OriginDH_y;
    Point1 = Originy*OriginRot_Negz*Point1DH;
    Point2 = Point1*Point2DH;
    Point3 = Point2*Point3DH;
    Point4 = Point3*Point4DH;

    % Extract Position in World Frame
    WorldOrigin = Originy(1:3,4);
    Point1Pos = Point1(1:3,4);
    Point2Pos = Point2(1:3,4);
    Point3Pos = Point3(1:3,4);
    Point4Pos = Point4(1:3,4);
   
    
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
    addOrientationArrows(Originx,scale);
    addOrientationArrows(Originy,scale);
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

