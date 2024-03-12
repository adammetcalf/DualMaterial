close all;
clear;
clc;

% oppositely magnetised
%load workspaceVariables.mat

% Dual material
load workspaceVariables2.mat

% Get S2
PointsmmS2(:,1) = 0 - PointsmmS(:,1);
PointsmmS2(:,3) = PointsmmS(:,3);

% Get J2
PointsmmJ2(:,1) = 0 - PointsmmJ(:,1);
PointsmmJ2(:,3) = PointsmmJ(:,3);

% Get C2
PointsmmC2(:,1) = 0 - PointsmmC(:,1);
PointsmmC2(:,3) = PointsmmC(:,3);


BaseTip = getTipVector(PointsmmB);
JTip = getTipVector(PointsmmJ);
J2Tip = getTipVector(PointsmmJ2);
CTip = getTipVector(PointsmmC);
STip = getTipVector(PointsmmS);

% Extract the tip positions as the final values in each array
BasePos = PointsmmB(end, :);
JPos = PointsmmJ(end, :);
J2Pos = PointsmmJ2(end, :);
CPos = PointsmmC(end, :);
SPos = PointsmmS(end, :);

% create Tip Iorientation and positions Base through J --- Dual Material
[VecArrayJ, PositionsJ] = interpolateJthroughJ(PointsmmB,PointsmmJ,PointsmmJ2);

% create Tip Iorientation and positions Base through S -- Dual Material
[VecArrayS, PositionsS] = interpolateSthroughS(PointsmmB,PointsmmS2,PointsmmS);

% create Tip Iorientation and positions Base J through S  --- Opp
% magentised
[VecArray, Positions] = interpolateJthroughS(PointsmmB,PointsmmJ,PointsmmS);

% create Tip Iorientation and positions J through C
[VecArrayJC, PositionsJC] = interpolateJthroughC(PointsmmJ, PointsmmC);

% create Tip Iorientation and positions S through C
[VecArraySC, PositionsSC] = interpolateJthroughC(PointsmmS, PointsmmC);

%% Dual Material Vector field

% For the XZ plane projection, we use the x and z components of these positions
%x = [PositionsJ(:,1);PositionsS(:,1);PositionsJC(:,1);PositionsSC(:,1)]; % x positions
%z = [PositionsJ(:,3);PositionsS(:,3);PositionsJC(:,3);PositionsSC(:,3)]; % z positions

% And the orientation unit vectors for the XZ plane
%u = [VecArrayJ(:,1);VecArrayS(:,1);VecArrayJC(:,1);VecArraySC(:,1)]; % x components of unit vectors
%w = [VecArrayJ(:,3);VecArrayS(:,3);VecArrayJC(:,3);VecArraySC(:,3)]; % z components of unit vectors


%% Oppositely Magnetised Vector Field

% For the XZ plane projection, we use the x and z components of these positions
x = [Positions(:,1);PositionsJC(:,1);PositionsSC(:,1)]; % x positions
z = [Positions(:,3);PositionsJC(:,3);PositionsSC(:,3)]; % z positions

% And the orientation unit vectors for the XZ plane
u = [VecArray(:,1);VecArrayJC(:,1);VecArraySC(:,1)]; % x components of unit vectors
w = [VecArray(:,3);VecArrayJC(:,3);VecArraySC(:,3)]; % z components of unit vectors


% Plot the vector field
figure(1);
quiver(x, z, u, w, 'AutoScale', 'off');
axis equal;
title('Tip Orientation Vector Field (XZ Plane)');
xlabel('X Position');
ylabel('Z Position'); % Note: This is plotted on the figure's y-axis for visualization
xlim([-25 25]);
ylim([-50 0]);

[Angles,Positions] = GetOrientationFromVertical(u,w,x,z);

% plot the tip orientation heatmap
figure(2)
scatter(Positions(:,1), Positions(:,3), 40, Angles, 'filled'); % Adjust the size (40) as needed

% Colorbar and colormap
colormap(jet); % You can choose any colormap you like (e.g., 'hsv' for a cyclic color scheme)
colorbar;
caxis([min(Angles), max(Angles)]); % Adjust the color limits to the range of your angles

% Labeling
xlabel('X Position');
ylabel('Y Position');
title('Position-Angle Plot');


%% Functions

% Get Unit Vector for tip orientation
function vec = getTipVector(Array)

    % flip array vertically
    Array = flip(Array,1);

    Tip = Array(1,:);
    PreTip = Array(2,:);

    vec  = Tip - PreTip;
    magnitude = sqrt(sum(vec.^2));

    vec = vec/magnitude;
end

% Sweep J to J2 tip through Base
function [VecArray, Positions] = interpolateJthroughJ(PointsmmB,PointsmmJ,PointsmmJ2)

    %% Flip arrays to get bottoms to top

    % flip array vertically
    PointsmmB = flip(PointsmmB,1);

    % flip array vertically
    PointsmmJ = flip(PointsmmJ,1);

    % flip array vertically
    PointsmmJ2 = flip(PointsmmJ2,1);

    %% Quadratic interpolate between tip positions

    % Combine the x and y tip coordinates into vectors
    x = [PointsmmJ(1,1), PointsmmB(1,1), PointsmmJ2(1,1)];
    y = [PointsmmJ(1,3), PointsmmB(1,3), PointsmmJ2(1,3)];
    
    % Fit a quadratic curve to the points
    p = polyfit(x, y, 2); 
    
    % Generate points along the fitted curve for plotting
    xFitTip = linspace(min(x), max(x), 101); % Generate 101 points for a smooth curve with defined centrepoint
    yFitTip = polyval(p, xFitTip); % Evaluate the polynomial at the xFit points

    %Reseize to encompass only between J and Centre
    xFitTip = xFitTip(1:52);
    yFitTip = yFitTip(1:52);

    %% Quadratic interpolate between Pretip positions

    % Combine the x and y tip coordinates into vectors
    x = [PointsmmJ(2,1), PointsmmB(2,1), PointsmmJ2(2,1)];
    y = [PointsmmJ(2,3), PointsmmB(2,3), PointsmmJ2(2,3)];
    
    % Fit a quadratic curve to the points
    p = polyfit(x, y, 2); 
    
    % Generate points along the fitted curve for plotting
    xFitPreTip = linspace(min(x), max(x), 101); % Generate 101 points for a smooth curve with defined centrepoint
    yFitPreTip = polyval(p, xFitPreTip); % Evaluate the polynomial at the xFit points

    %Reseize to encompass only between J and Centre
    xFitPreTip = xFitPreTip(1:52);
    yFitPreTip = yFitPreTip(1:52);

    %% Create vector of Unit vectors and vector of corresponding tip positions
    for i =1:length(xFitTip)
        
        Array = [xFitPreTip(i),0,yFitPreTip(i);xFitTip(i),0,yFitTip(i)];
        VecArray(i,:) =  getTipVector(Array);
        Positions(i,:) = [xFitTip(i),0,yFitTip(i)];
    end

end

% Sweep S to S2 tip through Base
function [VecArray, Positions] = interpolateSthroughS(PointsmmB,PointsmmJ,PointsmmJ2)

    %% Flip arrays to get bottoms to top

    % flip array vertically
    PointsmmB = flip(PointsmmB,1);

    % flip array vertically
    PointsmmJ = flip(PointsmmJ,1);

    % flip array vertically
    PointsmmJ2 = flip(PointsmmJ2,1);

    %% Quadratic interpolate between tip positions

    % Combine the x and y tip coordinates into vectors
    x = [PointsmmJ(1,1), PointsmmB(1,1), PointsmmJ2(1,1)];
    y = [PointsmmJ(1,3), PointsmmB(1,3), PointsmmJ2(1,3)];
    
    % Fit a quadratic curve to the points
    p = polyfit(x, y, 2); 
    
    % Generate points along the fitted curve for plotting
    xFitTip = linspace(min(x), max(x), 71); % Generate 101 points for a smooth curve with defined centrepoint
    yFitTip = polyval(p, xFitTip); % Evaluate the polynomial at the xFit points

    %Reseize to encompass only between J and Centre
    xFitTip = xFitTip(36:end);
    yFitTip = yFitTip(36:end);

    %% Quadratic interpolate between Pretip positions

    % Combine the x and y tip coordinates into vectors
    x2 = [PointsmmJ(2,1), PointsmmB(2,1), PointsmmJ2(2,1)];
    y2 = [PointsmmJ(2,3), PointsmmB(2,3), PointsmmJ2(2,3)];
    
    % Fit a quadratic curve to the points
    p = polyfit(x2, y2, 2); 
    
    % Generate points along the fitted curve for plotting
    xFitPreTip = linspace(min(x2), max(x2), 71); % Generate 101 points for a smooth curve with defined centrepoint
    yFitPreTip = polyval(p, xFitPreTip); % Evaluate the polynomial at the xFit points

    %Reseize to encompass only between J and Centre
    xFitPreTip = xFitPreTip(36:end);
    yFitPreTip = yFitPreTip(36:end);

    %% Create vector of Unit vectors and vector of corresponding tip positions
    for i =1:length(xFitTip)
        
        Array = [xFitPreTip(i),0,yFitPreTip(i);xFitTip(i),0,yFitTip(i)];
        VecArray(i,:) =  getTipVector(Array);
        Positions(i,:) = [xFitTip(i),0,yFitTip(i)];
    end

end

% Linearly interpolate between J tip and C tip
function [VecArray, Positions] = interpolateJthroughC(PointsmmJ, PointsmmC)

    %% Flip arrays to get bottoms to top

    % flip array vertically
    PointsmmJ = flip(PointsmmJ,1);

    % flip array vertically
    PointsmmC = flip(PointsmmC,1);

    %% linearly interpolate between tip positions
    
    % Generate points along the fitted curve for plotting
    xFitTip = linspace(PointsmmJ(1,1), PointsmmC(1,1), 30); % Generate 101 points for a smooth curve with defined centrepoint
    yFitTip = linspace(PointsmmJ(1,3), PointsmmC(1,3), 30); % Generate 101 points for a smooth curve with defined centrepoint

    %% linearly interpolate between tip positions
    
    % Generate points along the fitted curve for plotting
    xFitPreTip = linspace(PointsmmJ(2,1), PointsmmC(2,1), 30); % Generate 101 points for a smooth curve with defined centrepoint
    yFitPreTip = linspace(PointsmmJ(2,3), PointsmmC(2,3), 30); % Generate 101 points for a smooth curve with defined centrepoint

    %% Create vector of Unit vectors and vector of corresponding tip positions
    for i =1:length(xFitTip)
        
        Array = [xFitPreTip(i),0,yFitPreTip(i);xFitTip(i),0,yFitTip(i)];
        VecArray(i,:) =  getTipVector(Array);
        Positions(i,:) = [xFitTip(i),0,yFitTip(i)];
    end

end

% Sweep S to S2 tip through Base
function [VecArray, Positions] = interpolateJthroughS(PointsmmB,PointsmmJ,PointsmmS)

    %% Flip arrays to get bottoms to top

    % flip array vertically
    PointsmmB = flip(PointsmmB,1);

    % flip array vertically
    PointsmmJ = flip(PointsmmJ,1);

    % flip array vertically
    PointsmmS = flip(PointsmmS,1);

    %% Quadratic interpolate between tip positions

    % Combine the x and y tip coordinates into vectors
    x = [PointsmmJ(1,1), PointsmmB(1,1), PointsmmS(1,1)];
    y = [PointsmmJ(1,3), PointsmmB(1,3), PointsmmS(1,3)];
    
    % Fit a quadratic curve to the points
    p = polyfit(x, y, 2); 
    
    % Generate points along the fitted curve for plotting
    xFitTip = linspace(min(x), max(x), 71); % Generate 101 points for a smooth curve with defined centrepoint
    yFitTip = polyval(p, xFitTip); % Evaluate the polynomial at the xFit points

    %Reseize to encompass only between J and Centre
    %xFitTip = xFitTip(36:end);
    %yFitTip = yFitTip(36:end);

    %% Quadratic interpolate between Pretip positions

    % Combine the x and y tip coordinates into vectors
    x2 = [PointsmmJ(2,1), PointsmmB(2,1), PointsmmS(2,1)];
    y2 = [PointsmmJ(2,3), PointsmmB(2,3), PointsmmS(2,3)];
    
    % Fit a quadratic curve to the points
    p = polyfit(x2, y2, 2); 
    
    % Generate points along the fitted curve for plotting
    xFitPreTip = linspace(min(x2), max(x2), 71); % Generate 101 points for a smooth curve with defined centrepoint
    yFitPreTip = polyval(p, xFitPreTip); % Evaluate the polynomial at the xFit points

    %Reseize to encompass only between J and Centre
    %xFitPreTip = xFitPreTip(36:end);
    %yFitPreTip = yFitPreTip(36:end);

    %% Create vector of Unit vectors and vector of corresponding tip positions
    for i =1:length(xFitTip)
        
        Array = [xFitPreTip(i),0,yFitPreTip(i);xFitTip(i),0,yFitTip(i)];
        VecArray(i,:) =  getTipVector(Array);
        Positions(i,:) = [xFitTip(i),0,yFitTip(i)];
    end

end

function [Angles,Positions] = GetOrientationFromVertical(u,w,x,z)

    % Evaluate Vector angles from negative z
    angleOrigin = [0,0,-1];

    % Calculate the magnitude angleOrigin
    normangleOrigin = norm(angleOrigin);

    % recombine Vectors
    for i = 1:length(u)
        vectors(i,1) = u(i);
        vectors(i,2) = 0;
        vectors(i,3) = w(i);
    
        Positions(i,1) = x(i);
        Positions(i,2) = 0;
        Positions(i,3) = z(i);

        % Calculate the dot product of vec and vecB
        dotProduct = dot(vectors(i,:), angleOrigin);

        % Calculate the magnitude (norms) of vec and vecB
        normVec = norm(vectors(i,:));

        % Calculate the cosine of the angle
        cosTheta = dotProduct / (normVec * normangleOrigin);
        
        % Calculate the angle in radians
        angleRadians = acos(cosTheta);

        Angles(i) = rad2deg(angleRadians);
    end

end