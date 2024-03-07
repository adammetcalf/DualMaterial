close all;
clear;
clc;

% Rod is composed of N+1 rigid links joined by N flexible joints. Each
% joint has 2Dof for bending (ignoring extensions, shear and twist of the
% rod) The deformation of a signle rod may be described by the rotation
% angles Theta1(k) and Theta2(k). For this 2D analysis only Theta2 matters.

%% Define Tentacle Major components

% Tentacle Length (40mm)
length = 40/1000; %(m)

%Number of links
N = 10;

%Link Length
L = length/N; %(m)

% Link Radii (m)
radius = 1e-03;

% EcoFlex 0030 Density (kg m^-3)
rho = 1070; 

% Youngs Modulus Ecoflex 0030 (125 kPa)
E = 125000;

% Initialise Material properties (for each segment)
[Mass,I,k,~] = MaterialProperties(radius,rho,L,E);

%% Initialise angles and forces

% initialise Theta2 angles
for i=1:N
    T2(i,1) = 0;
end

% initialise forces
for i=1:N
    Forcex(1,i) = 0;
    Forcey(1,i) = 0;
    Forcez(1,i) = 0;
    TorqueX(1,i) = 0;
    TorqueY(1,i) = 0;
    TorqueZ(1,i) = 0;
end
Force = [Forcex; Forcey; Forcez; TorqueX; TorqueY; TorqueZ];

% F=kx, Tau =ktheta. Not that here we care about Fx, Fz and TauY


%% Define Tentacle distal end magnetic moments.


p = forwardKinematics(T2,L);

plot3([0,p(1,1:5)],[0,p(2,1:5)],[0,p(3,1:5)],...
    'bo-', 'MarkerSize', 4, 'MarkerFaceColor', 'b','LineWidth', 1);
hold on
plot3([p(1,5:6)],[p(2,5:6)],[p(3,5:6)],...
    'r-','MarkerFaceColor', 'r','LineWidth', 1);
plot3([p(1,6:end)],[p(2,6:end)],[p(3,6:end)],...
    'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1);
axis([-0.05 0.05 -0.05 0.05 -0.05 0.05]);
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')
grid on

%% Functions

% forward kinematics
function p = forwardKinematics(T2,L)

    for i = 1:length(T2)
        if i==1
            p(:,i) = [cos(T2(i))*cos(0)*L
            cos(T2(i))*sin(0)*L
            sin(T2(i))*L];
            
        else
           p(:,i) = p(:,i-1) + [cos(T2(i))*cos(0)*L
           cos(T2(i))*sin(0)*L
           sin(T2(i))*L];
        end
    end

end

% Adjust Angle
function T2 = AdjustAngle(T2,index,Angle)

Angle = deg2rad(Angle);

    for i=1:length(T2)
        if i >= index
            T2(i) = T2(i) + Angle;
        end
    end
end