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

theta0 = deg2rad(0);
alpha0 = deg2rad(90); %% Initial starting position defined by this value.

OriginDH = [cos(theta0) -sin(theta0)*cos(alpha0) sin(theta0)*sin(alpha0) 0
            sin(theta0) cos(theta0)*cos(alpha0) -cos(theta0)*sin(alpha0) 0
            0 sin(alpha0) cos(alpha0) 0
            0 0 0 1];

theta1 = deg2rad(0);
alpha1 = deg2rad(0);

Point1DH = [cos(theta1) -sin(theta1)*cos(alpha1) sin(theta1)*sin(alpha1) 0
            sin(theta1) cos(theta1)*cos(alpha1) -cos(theta1)*sin(alpha1) 0
            0 sin(alpha1) cos(alpha1) 0.01
            0 0 0 1];

theta2 = deg2rad(0);
alpha2 = deg2rad(0);

Point2DH = [cos(theta2) -sin(theta2)*cos(alpha2) sin(theta2)*sin(alpha2) 0
            sin(theta2) cos(theta2)*cos(alpha2) -cos(theta2)*sin(alpha2) 0
            0 sin(alpha2) cos(alpha2) 0.01
            0 0 0 1];

theta3 = deg2rad(0);
alpha3 = deg2rad(0);

Point3DH = [cos(theta3) -sin(theta3)*cos(alpha3) sin(theta3)*sin(alpha3) 0
            sin(theta3) cos(theta3)*cos(alpha3) -cos(theta3)*sin(alpha3) 0
            0 sin(alpha3) cos(alpha3) 0.01
            0 0 0 1];

theta4 = deg2rad(0);
alpha4 = deg2rad(0);

Point4DH = [cos(theta4) -sin(theta4)*cos(alpha4) sin(theta4)*sin(alpha4) 0
            sin(theta4) cos(theta4)*cos(alpha4) -cos(theta4)*sin(alpha4) 0
            0 sin(alpha4) cos(alpha4) 0.01
            0 0 0 1];


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

% Link Lengths (m)
Length01 = Point1DH(3,4);
Length12 = Point2DH(3,4);
Length23 = Point3DH(3,4);
Length34 = Point4DH(3,4);

% Link Volumes (m^3)
Vol01 = Length01*pi*(radius^2);
Vol12 = Length12*pi*(radius^2);
Vol23 = Length23*pi*(radius^2);
Vol34 = Length34*pi*(radius^2);

% Link Masses (consider as a point mass, with the mass of each link acting
% at the endpoint - TODO: Adjust this to act at the COM of each link.)
Mass1 = rho*Vol01;
Mass2 = rho*Vol12;
Mass3 = rho*Vol23;
Mass4 = rho*Vol34;

% Define Inertia (kg m^-2) 
%%% TO DO: Inertia evaluated about correct position?
%%% SOLID CYCLNDER ABOUT ENDPOINT DIAMETER
%%% I = 1/4*m*R^2 + 1/3*m*L^2

I1 = (1/4)*Mass1*radius^2 + (1/3)*Mass1*Length01;
I2 = (1/4)*Mass2*radius^2 + (1/3)*Mass2*Length12;
I3 = (1/4)*Mass3*radius^2 + (1/3)*Mass3*Length23;
I4 = (1/4)*Mass4*radius^2 + (1/3)*Mass4*Length34;

% Define Joint Stiffness (Nm rad^-1)
% k_m = EI L^-1

% NOTE, stiffness evaluated for Link01 is applied at joint0 etc

k_1 = (E*I1)/Length01;
k_2 = (E*I2)/Length12;
k_3 = (E*I3)/Length23;
k_4 = (E*I4)/Length34;


%% Magnetic Properties of the tentacle

% Define Magnetic Moments (Magnitudes)
m1 = 0.05;                             % magnitude of magnetic moment of point 1
m2 = 0.05;                             % magnitude of magnetic moment of point 2
m3 = 0.05;                             % magnitude of magnetic moment of point 3
m4 = 0.05;                             % magnitude of magnetic moment of point 4

%% Place Tentacle
Origin = StartTransform*OriginDH;
Point1 = StartTransform*OriginDH*Point1DH;
Point2 = StartTransform*OriginDH*Point1DH*Point2DH;
Point3 = StartTransform*OriginDH*Point1DH*Point2DH*Point3DH;
Point4 = StartTransform*OriginDH*Point1DH*Point2DH*Point3DH*Point4DH;

WorldOrigin = Origin(1:3,4);
Point1Pos = Point1(1:3,4);
Point2Pos = Point2(1:3,4);
Point3Pos = Point3(1:3,4);
Point4Pos = Point4(1:3,4);

%% Plot Tentacle
figure(1)
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