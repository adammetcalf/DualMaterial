close all;
clear;
clc;

% Pure NdFEB Script

%% Define the world

WorldOrigin = [0,0,0];
Gravity = [0,0,-9.81];                  % m/s^2

%% Define the rigid body points

Point1 = [0,-0.01,0];                   % initial condition for point 1
Point2 = [0,-0.02,0];                   % initial condition for point 2
Point3 = [0,-0.03,0];                   % initial condition for point 3
Point4 = [0,-0.04,0];                   % initial condition for point 3
Mass = 0.001;                           % kg
I = 0.002;                              % kg/m^2

m1 = 0.05;                             % magnitude of magnetic moment of point 1
m2 = 0.05;                             % magnitude of magnetic moment of point 2
m3 = 0.05;                             % magnitude of magnetic moment of point 3
m4 = 0.05;                             % magnitude of magnetic moment of point 3

% Initialise Velocities                 % m/s
Velocity1 = [0,0,0];
Velocity2 = [0,0,0];
Velocity3 = [0,0,0];
Velocity4 = [0,0,0];

% Initialize angular velocities         %Rad/s
Omega1 = [0, 0, 0];
Omega2 = [0, 0, 0];
Omega3 = [0, 0, 0];
Omega4 = [0, 0, 0];


%% Define Damping Factor
DampingFactor = 0.02; 
displacementFactor = 0.005;

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
totalTime = 1; 
Iterations = totalTime/dt;

%% Define Magnetic Stuff
mu0 = 4*pi*1e-7; % Permeability of free space
B = [0,0,25e-3]; % Magnetic field strength in Tesla (25 mT) in +z direction

% Create the magnetic field
xLow = -0.05;   yLow = -0.05;   zLow = -0.1;
xHigh = 0.05;   yHigh = 0.05;   zHigh = 0.1;
xIncr = 0.005;  yIncr = 0.005;  zIncr = 0.005; %steps size

[x, y, z] = meshgrid(xLow:xIncr:xHigh, yLow:yIncr:yHigh, zLow:zIncr:zHigh);

BxCoil = B(1)*ones(size(x));
ByCoil = B(2)*ones(size(y));
BzCoil = B(3)*ones(size(z));

threshold = 0.001; %Singularity threshold

%% Simulation Loop
for step = 1:Iterations

    disp(step)

    m1_dir = getDirection(WorldOrigin, Point1);
    m2_dir = getDirection(Point1, Point2);
    m3_dir = getDirection(Point2, Point3);
    m4_dir = getDirection(Point3, Point4);

    % Rotate magnetic moment directions
    m1_dir = rotateVector(m1_dir, Omega1, norm(Omega1) * dt);
    m2_dir = rotateVector(m2_dir, Omega2, norm(Omega2) * dt);
    m3_dir = rotateVector(m3_dir, Omega3, norm(Omega3) * dt);
    m4_dir = rotateVector(m4_dir, Omega4, norm(Omega4) * dt);

    % Calculate Magnetic Moments
    m1_Vector = m1 * m1_dir; % m1 magnitude along Link1's direction
    m2_Vector = m2 * m2_dir; % m2 magnitude along Link1's direction
    m3_Vector = m3 * m3_dir; % m3 magnitude along Link1's direction
    m4_Vector = m4 * m4_dir; % m4 magnitude along Link1's direction

    % Calculate field components for the dipoles
    [Bx1,By1,Bz1] = evaluateField(x,y,z,Point1,mu0,m1_Vector,threshold);
    [Bx2,By2,Bz2] = evaluateField(x,y,z,Point2,mu0,m2_Vector,threshold);
    [Bx3,By3,Bz3] = evaluateField(x,y,z,Point3,mu0,m3_Vector,threshold);
    [Bx4,By4,Bz4] = evaluateField(x,y,z,Point4,mu0,m4_Vector,threshold);

    % Sum the magnetic fields from both dipoles
    Bx_total = BxCoil + Bx1 + Bx2 + Bx3 + Bx4;
    By_total = ByCoil + By1 + By2 + By3 + By4;
    Bz_total = BzCoil + Bz1 + Bz2 + Bz3 + Bz4;

    % Find meshgrid indices for each dipole
    [idx1, idy1, idz1] = findClosestGridPoint(x, y, z, Point1);
    [idx2, idy2, idz2] = findClosestGridPoint(x, y, z, Point2);
    [idx3, idy3, idz3] = findClosestGridPoint(x, y, z, Point3);
    [idx4, idy4, idz4] = findClosestGridPoint(x, y, z, Point4);
    
    % Calculate magnetic forces
    F1 = f_getForce(Bx_total, By_total, Bz_total, idx1, idy1, idz1, m1_Vector, xIncr, yIncr, zIncr);
    F2 = f_getForce(Bx_total, By_total, Bz_total, idx2, idy2, idz2, m2_Vector, xIncr, yIncr, zIncr);
    F3 = f_getForce(Bx_total, By_total, Bz_total, idx3, idy3, idz3, m3_Vector, xIncr, yIncr, zIncr);
    F4 = f_getForce(Bx_total, By_total, Bz_total, idx4, idy4, idz4, m4_Vector, xIncr, yIncr, zIncr);

    % Calculate magnetic torques
    % Calculate Torques Point1
    T1_Coil = f_getTorque(BxCoil,ByCoil,BzCoil,idx1, idy1, idz1, m1_Vector);
    T1_2 = f_getTorque(Bx2,By2,Bz2,idx1, idy1, idz1, m1_Vector);
    T1_3 = f_getTorque(Bx3,By3,Bz3,idx1, idy1, idz1, m1_Vector);
    T1_4 = f_getTorque(Bx4,By4,Bz4,idx1, idy1, idz1, m1_Vector);

    % Calculate Torques Point2
    T2_Coil = f_getTorque(BxCoil,ByCoil,BzCoil,idx2, idy2, idz2, m2_Vector);
    T2_1 = f_getTorque(Bx1,By1,Bz1,idx2, idy2, idz2, m2_Vector);
    T2_3 = f_getTorque(Bx3,By3,Bz3,idx2, idy2, idz2, m2_Vector);
    T2_4 = f_getTorque(Bx4,By4,Bz4,idx2, idy2, idz2, m2_Vector);

    % Calculate Torques Point3
    T3_Coil = f_getTorque(BxCoil,ByCoil,BzCoil,idx3, idy3, idz3, m3_Vector);
    T3_1 = f_getTorque(Bx1,By1,Bz1,idx3, idy3, idz3, m3_Vector);
    T3_2 = f_getTorque(Bx2,By2,Bz2,idx3, idy3, idz3, m3_Vector);
    T3_4 = f_getTorque(Bx4,By4,Bz4,idx3, idy3, idz3, m3_Vector);    

    % Calculate Torques Point4
    T4_Coil = f_getTorque(BxCoil,ByCoil,BzCoil,idx4, idy4, idz4, m4_Vector);
    T4_1 = f_getTorque(Bx1,By1,Bz1,idx4, idy4, idz4, m4_Vector);
    T4_2 = f_getTorque(Bx2,By2,Bz2,idx4, idy4, idz4, m4_Vector);
    T4_3 = f_getTorque(Bx3,By3,Bz3,idx4, idy4, idz4, m4_Vector); 

    % Sum the torques to get the total torque on dipoles
    torque1 = T1_Coil+ T1_2 + T1_3+ T1_4; % torque acting on Point1
    torque2 = T2_Coil+ T2_1 + T2_3+ T2_4; % torque acting on Point2
    torque3 = T3_Coil+ T3_1 + T3_2+ T3_4; % torque acting on Point3
    torque4 = T4_Coil+ T4_1 + T4_2+ T4_3; % torque acting on Point4

    clear T1_Coil T2_Coil T3_Coil T4_Coil T1_2 T1_3 T1_4 T2_1 T2_3 T2_4 T3_1 T3_2 T3_4 T4_1  T4_2 T4_3
    
    % Calculate accelerations due to magnetic forces
    acceleration1 = F1 / Mass + Gravity;
    acceleration2 = F2 / Mass + Gravity;
    acceleration3 = F3 / Mass + Gravity;
    acceleration4 = F4 / Mass + Gravity;

    Velocity1 = Velocity1 + (acceleration1 - DampingFactor * Velocity1) * dt;
    Velocity2 = Velocity2 + (acceleration1 - DampingFactor * Velocity2) * dt;
    Velocity3 = Velocity3 + (acceleration1 - DampingFactor * Velocity3) * dt;
    Velocity4 = Velocity4 + (acceleration1 - DampingFactor * Velocity4) * dt;

    %% Update angular velocities based on torques
    Omega1 = Omega1 + (torque1 / I) * dt;
    Omega2 = Omega2 + (torque2 / I) * dt;
    Omega3 = Omega3 + (torque3 / I) * dt;
    Omega4 = Omega4 + (torque4 / I) * dt;

    for i = 1:numel(Point1)
        if norm(Omega1) > 0  % Assuming movement is influenced by the rate of rotation
            Point1(i) = Point1(i) + sign(Omega1(i)) * displacementFactor; % displacementFactor is a conceptual scaling factor
        end
    end

    for i = 1:numel(Point2)
        if norm(Omega2) > 0  % Assuming movement is influenced by the rate of rotation
            Point2(i) = Point2(i) + sign(Omega2(i)) * displacementFactor; % displacementFactor is a conceptual scaling factor
        end
    end

    for i = 1:numel(Point3)
        if norm(Omega3) > 0  % Assuming movement is influenced by the rate of rotation
            Point3(i) = Point3(i) + sign(Omega3(i)) * displacementFactor; % displacementFactor is a conceptual scaling factor
        end
    end

    for i = 1:numel(Point4)
        if norm(Omega4) > 0  % Assuming movement is influenced by the rate of rotation
            Point4(i) = Point4(i) + sign(Omega4(i)) * displacementFactor; % displacementFactor is a conceptual scaling factor
        end
    end    


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

    % Plot the current state
    clf; % Clear the current figure
    plot3([WorldOrigin(1), Point1(1), Point2(1), Point3(1), Point4(1)],...
          [WorldOrigin(2), Point1(2), Point2(2), Point3(2), Point4(2)],...
          [WorldOrigin(3), Point1(3), Point2(3), Point3(3), Point4(3)], 'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1);
    hold on
    %plot Magnetic field
    quiver3(x, y, z, Bx_total, By_total, Bz_total, 0.75, 'b','LineWidth', 0.5);
    hold off;
    axis equal;
    axis([-0.05 0.05 -0.05 0.05 -0.1 0.1]);
    % Set the figure to full screen using WindowState
    set(gcf, 'WindowState', 'maximized');
    drawnow;
end

%% Functions
function Pos2 = enforceLinkConstraint(P1, P2, L2)
    % Enforce the distance between Point1 and Point2
    dir = (P2 - P1) / norm(P2 - P1); % Direction from Point1 to Point2
    Pos2 = P1 + dir * L2; % Adjust Point2 position correctly
end

function dir = getDirection(P1, P2)
    dir = (P2 - P1) / norm(P2 - P1); % Direction from Point1 to Point2;
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
    dBy_dy = (By(x_index, y_index_plus, z_index) - By(x_index, y_index_minus, z_index)) / (2*dy);
    dBz_dz = (Bz(x_index, y_index, z_index_plus) - Bz(x_index, y_index, z_index_minus)) / (2*dz);

    % Calculate the force
    force = [m(1) * dBx_dx, m(2) * dBy_dy, m(3) * dBz_dz];

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

% Function to find the closest grid point
function [idx, idy, idz] = findClosestGridPoint(x, y, z, p)
    [~, idx] = min(abs(x(1,:,1) - p(1)));
    [~, idy] = min(abs(y(:,1,1) - p(2)));
    [~, idz] = min(abs(z(1,1,:) - p(3)));
end

function v_rot = rotateVector(v, k, theta)
    % Rotate vector v around unit vector k by angle theta (in radians)
    % Using Rodrigues' rotation formula

    % Check if the magnitude of the angular velocity vector is zero
    if norm(k) == 0
        % If yes, return the original vector without any changes
        v_rot = v;
    else
        % Otherwise, proceed with the rotation
        k = k / norm(k); % Ensure k is a unit vector
        v_rot = v*cos(theta) + cross(k, v)*sin(theta) + k*dot(k, v)*(1 - cos(theta));
    end
end