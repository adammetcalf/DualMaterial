close all;
clear;
clc;

%Load images
imgBase = imread('C:\Users\bsgx043\Desktop\FEBaseNdDistal\BaseThresh.png');
imgS = imread('C:\Users\bsgx043\Desktop\FEBaseNdDistal\3. MagR\SShapeThresh.png');
imgC = imread('C:\Users\bsgx043\Desktop\FEBaseNdDistal\2. FieldRtoL\CShapeThresh.png');
imgJ = imread('C:\Users\bsgx043\Desktop\FEBaseNdDistal\2. FieldRtoL\JShapeThresh.png');

% Define number of points to allocate
points = 50;

%% Display
%figure(1) 
%imshow(imgBase)

%figure(2) 
%imshow(imgS)

%figure(3) 
%imshow(imgC)

%figure(4) 
%imshow(imgJ)

%% Use Base image to get mm of each pixel

% Define centreline
[x,y] = getCentreline(imgBase);
[x,y] = extendCentreLineUp(x,y,imgBase);
[x,y] = extendCentreLineDown(x,y,imgBase);

% Get equidistant locations
XYLength = length(x);

% Generate x equidistant indices including the first and last
indices = round(linspace(1, XYLength, points));

% Select elements from A using the indices
XPoints = x(indices);
YPoints = y(indices);

% 40mm length of tentacle.
VecLength = [x(end) - x(1), y(end)-y(1)];
VecLength = sqrt(VecLength(1)^2 + VecLength(2)^2); %number of pixels per 40mm
VecLength = VecLength/40; % Number of pixels/mm
VecLength = 1/VecLength; % mm/pixel


%% Analyse Images
img = imgBase;
[xB,zB,XPointsB,ZPointsB,PointsmmB,CentremmB,sweepXB,sweepYB,sweepZB] = analyseImage(img,VecLength,points);
img = imgJ;
[xJ,zJ,XPointsJ,ZPointsJ,PointsmmJ,CentremmJ,sweepXJ,sweepYJ,sweepZJ] = analyseImage(img,VecLength,points);
img = flip(imgJ,2);
[xJ2,zJ2,XPointsJ2,ZPointsJ2,PointsmmJ2,CentremmJ2,sweepXJ2,sweepYJ2,sweepZJ2] = analyseImage(img,VecLength,points);
img = imgS;
[xS,zS,XPointsS,ZPointsS,PointsmmS,CentremmS,sweepXS,sweepYS,sweepZS] = analyseImage(img,VecLength,points);
img = flip(imgS,2);
%[xS2,zS2,XPointsS2,ZPointsS2,PointsmmS2,CentremmS2,sweepXS2,sweepYS2,sweepZS2] = analyseImage(img,VecLength,points);
PointsmmS2(:,1) = 0 - PointsmmS(:,1);
PointsmmS2(:,3) = PointsmmS(:,3);
CentremmS2(:,1) = 0 - CentremmS(:,1);
CentremmS2(:,3) = CentremmS(:,3);
img = imgC;
[xC,zC,XPointsC,ZPointsC,PointsmmC,CentremmC,sweepXC,sweepYC,sweepZC] = analyseImage(img,VecLength,points);
img = flip(imgC,2);
%[xC2,zC2,XPointsC2,ZPointsC2,PointsmmC2,CentremmC2,sweepXC2,sweepYC2,sweepZC2] = analyseImage(img,VecLength,points);
PointsmmC2(:,1) = 0 - PointsmmC(:,1);
PointsmmC2(:,3) = PointsmmC(:,3);
CentremmC2(:,1) = 0 - CentremmC(:,1);
CentremmC2(:,3) = CentremmC(:,3);


%% Analyse Workspace
%PS = analyseWorkspace(CentremmB,CentremmJ2,CentremmJ); %works
%PS2 = analyseWorkspace(CentremmB,CentremmS2,CentremmS); %works
%PS3 = analyseWorkspaceC(CentremmB,CentremmC, CentremmS2,CentremmS);
[PS,PSx,PSz] = tipWorkspace(CentremmB,CentremmC,CentremmC2,CentremmS,CentremmS2,CentremmJ,CentremmJ2);
img = combineImages(imgBase,imgS,imgC,imgJ);
[X3D,Y3D,Z3D] = revolveTipWorkspace(PSx,PSz);


%% Plot
figure(5)
subplot(1,2,1)
imshow(img)
hold on
subplot(1,2,2)
imshow(img)
hold on
plot(XPointsB,ZPointsB,'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r')
plot(XPointsJ,ZPointsJ,'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r')
plot(XPointsC,ZPointsC,'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r')
plot(XPointsS,ZPointsS,'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r')



%% Plot
figure(6)
subplot(1,2,1)
plot(CentremmB(:,1),CentremmB(:,3),'k')
hold on
plot(CentremmJ(:,1),CentremmJ(:,3),'r')
plot(CentremmJ2(:,1),CentremmJ2(:,3),'r')
plot(CentremmC(:,1),CentremmC(:,3),'r')
plot(CentremmC2(:,1),CentremmC2(:,3),'r')
plot(CentremmS(:,1),CentremmS(:,3),'r')
%plot(CentremmS2(:,1),CentremmS2(:,3),'r-')
plot(PointsmmB(:,1),PointsmmB(:,3),'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'k')
plot(PointsmmJ(:,1),PointsmmJ(:,3),'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r')
plot(PointsmmJ2(:,1),PointsmmJ2(:,3),'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r')
plot(PointsmmC(:,1),PointsmmC(:,3),'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r')
plot(PointsmmC2(:,1),PointsmmC2(:,3),'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r')
plot(PointsmmS(:,1),PointsmmS(:,3),'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r')
plot(PointsmmS2(:,1),PointsmmS2(:,3),'ro-', 'MarkerSize', 2, 'MarkerFaceColor', 'r')
xlim([-25 25]);
ylim([-50 0]);
xlabel('X (mm)')
ylabel('Z (mm)')
axis equal
subplot(1,2,2)
plot(PointsmmB(:,1),PointsmmB(:,3),'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'k')
hold on
plot(PointsmmJ(:,1),PointsmmJ(:,3),'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r')
plot(PointsmmJ2(:,1),PointsmmJ2(:,3),'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r')
plot(PointsmmC(:,1),PointsmmC(:,3),'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r')
plot(PointsmmC2(:,1),PointsmmC2(:,3),'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r')
plot(PointsmmS(:,1),PointsmmS(:,3),'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r')
plot(PointsmmS2(:,1),PointsmmS2(:,3),'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r')
plot(PS)
xlim([-25 25]);
ylim([-50 0]);
xlabel('X (mm)')
ylabel('Z (mm)')
axis equal

figure(7)
% Plot the 3D volumetric shape
surf(X3D, Y3D, Z3D, 'FaceColor', 'red', 'EdgeColor', 'none','FaceAlpha', 0.5);
hold on
plot3(PointsmmB(:,1),PointsmmB(:,2),PointsmmB(:,3),'ko-', 'MarkerSize', 2, 'MarkerFaceColor', 'k')
plot3(PointsmmJ(:,1),PointsmmJ(:,2),PointsmmJ(:,3),'ko-', 'MarkerSize', 2, 'MarkerFaceColor', 'k')
plot3(PointsmmC(:,1),PointsmmC(:,2),PointsmmC(:,3),'ko-', 'MarkerSize', 2, 'MarkerFaceColor', 'k')
plot3(PointsmmS(:,1),PointsmmS(:,2),PointsmmS(:,3),'ko-', 'MarkerSize', 2, 'MarkerFaceColor', 'k')
axis equal; % Ensure aspect ratio is 1:1:1
xlabel('X mm');
ylabel('Y mm');
zlabel('Z mm');
xlim([-25 25]);
ylim([-25 25]);
zlim([-50 0]);
%camlight left; % Add some lighting for better visualization
%lighting phong; % Use Phong lighting for a smoother appearance


%% Functions

% Function to get the centreline
function [x,y] = getCentreline(img)

    % Convert image to binary if it's not already
    if size(img, 3) == 3
        imgBinary = imbinarize(rgb2gray(img));
    else
        imgBinary = imgBase > 0; % Assuming non-zero pixels are the rod
    end
    
    % Skeletonize the binary image to find the medial axis of the rod
    
    B = bwmorph(~imgBinary,'skel',inf); % skeletonize
    B = imfill(B,'holes'); % get rid of any tiny closed paths
    B = bwmorph(B,'skel',inf); % finish skeletonizing
    % remove endpoints until only two are left
    
    while true
        e = bwmorph(B,'endpoints');
        if nnz(e)>2
            B = B & ~e;
        else
            break;
        end
    end

    
    % Find the coordinates of the skeleton
    [y, x] = find(B);
    
    % Sort points by 'y' to ensure top-to-bottom ordering
    [~, order] = sort(y, 'ascend');
    x = x(order);
    y = y(order);

end

% Function to detect the edges
function [img,thinnedEdges] = detectEdges(img)

    % Convert the image to grayscale if it is not already
    if size(img, 3) == 3
        imgGray = rgb2gray(img);
    else
        imgGray = img;
    end
    
    % Step 1: Edge Detection using the Canny method
    edges = edge(imgGray, 'Canny');
    
    % Step 2: Thinning to get an approximation of the centerline
    thinnedEdges = bwmorph(edges, 'thin', Inf);
    
    
    % Overlay the thinned edges on the original image
    for i = 1:size(thinnedEdges, 1)
        for j = 1:size(thinnedEdges, 2)
            if thinnedEdges(i, j)
                img(i, j, :) = [0, 0, 255]; % Blue color for the centerline
            end
        end
    end

end


function [x,y] = extendCentreLineUp(x,y,img)
    
    % Convert the image to grayscale if it is not already
    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    %reverse Array directions
    x = flip(x);
    y = flip(y);
    

    % get vector between last 30 points (x2-x1,y2-y1) and normalise
    sizeXY = length(x);  
    %vec = [(x(sizeXY) - x(sizeXY-30))/30,(y(sizeXY) - y(sizeXY-30))/30];
    vec = [0,-1];

    stop = false;
    %define Increment
    inc = 1;

    while stop == false

        %place a new x and a new y at the end of their respective vectors,
        %using the deirection given by vec
        x(sizeXY+inc) = round(x(sizeXY)+inc*vec(1));
        y(sizeXY+inc) = round(y(sizeXY)+inc*vec(2));

        %check img to see if we have hit the white pixel yet
        if img(y(end),x(end)) == 255
            stop = true;
            y(end) = [];
            x(end) = [];
        else
            inc = inc+1;
        end
    end

    %return Array directions to origin inputs
    x = flip(x);
    y = flip(y);

end

function [x,y] = extendCentreLineDown(x,y,img)

    % Convert the image to grayscale if it is not already
    if size(img, 3) == 3
        img = rgb2gray(img);
    end
  
    PointLength = 10;
    % get vector between last 30 points (x2-x1,y2-y1) and normalise
    sizeXY = length(x);  
    vec = [(x(sizeXY) - x(sizeXY-PointLength))/PointLength,(y(sizeXY) - y(sizeXY-PointLength))/PointLength];


    stop = false;
    %define Increment
    inc = 1;

    while stop == false

        %place a new x and a new y at the end of their respective vectors,
        %using the deirection given by vec
        x(sizeXY+inc) = round(x(sizeXY)+inc*vec(1));
        y(sizeXY+inc) = round(y(sizeXY)+inc*vec(2));

        %check img to see if we have hit the white pixel yet
        if img(y(end),x(end)) == 255
            stop = true;
            y(end) = [];
            x(end) = [];
        else
            inc = inc+1;
        end
    end

end

function [x,z,XPoints,ZPoints,Pointsmm,Centremm,sweepX,sweepY,sweepZ] = analyseImage(img,VecLength,points)

    %% Skeleton method
    % Define centreline
    [x,z] = getCentreline(img);
    [x,z] = extendCentreLineUp(x,z,img);
    [x,z] = extendCentreLineDown(x,z,img);
    
    % Get equidistant locations
    XYLength = length(x);
    
    % Generate x equidistant indices including the first and last
    indices = round(linspace(1, XYLength, points));
    
    % Select elements from A using the indices
    XPoints = x(indices);
    ZPoints = z(indices);
    
    % Offset all points (in mm) where[x1,y1] is at [0mm,0mm]
    xmm = (x-x(1))*VecLength;
    zmm = -(z-z(1))*VecLength;
    ymm = zeros(length(xmm),1);
    XPointsmm = (XPoints-x(1))*VecLength;
    ZPointsmm = -(ZPoints-z(1))*VecLength;
    YPointsmm = zeros(length(XPointsmm),1);

    Centremm = [xmm,ymm,zmm];    
    Pointsmm = [XPointsmm,YPointsmm,ZPointsmm];
    
    % Number of points to generate around the Z-axis
    numPoints = 360;
    
    % Preallocate array for efficiency
    sweepData = zeros(size(Pointsmm,1) * numPoints, 3);
    
    % Sweep angle from 0 to 2*pi (360 degrees)
    theta = linspace(0, 2*pi, numPoints);
    
    % Loop through all angles to create the sweep
    for i = 1:numPoints
        Rz = [cos(theta(i)) -sin(theta(i)) 0; sin(theta(i)) cos(theta(i)) 0; 0 0 1]; % Rotation matrix around Z-axis
        rotatedPoints = (Rz * Pointsmm')'; % Apply rotation to all points
        sweepData((i-1)*size(Pointsmm,1)+1:i*size(Pointsmm,1), :) = rotatedPoints; % Store rotated points
    end
    
    % Reshape the data for plotting
    sweepX = reshape(sweepData(:,1), [size(Pointsmm,1), numPoints]);
    sweepY = reshape(sweepData(:,2), [size(Pointsmm,1), numPoints]);
    sweepZ = reshape(sweepData(:,3), [size(Pointsmm,1), numPoints]);  

end

function PS = analyseWorkspace(Centremm,Leftmm,Rightmm)

    TipEnd_left = Leftmm(end,:);
    TipEnd_centre = Centremm(end,:);
    TipEnd_right = Rightmm(end,:);
    
    
    % Combine the x and y coordinates into vectors
    x = [TipEnd_left(1), TipEnd_centre(1), TipEnd_right(1)];
    y = [TipEnd_left(3), TipEnd_centre(3), TipEnd_right(3)];
    
    % Fit a quadratic curve to the points
    p = polyfit(x, y, 2); 
    
    % Generate points along the fitted curve for plotting
    xFit = linspace(min(x), max(x), 100); % Generate 100 points for a smooth curve
    yFit = polyval(p, xFit); % Evaluate the polynomial at the xFit points
    
    % Create Geometric shape for workspace visualisation
    PSx = [Leftmm(:,1);xFit';flip(Rightmm(:,1))];
    PSy = [Leftmm(:,3);yFit';flip(Rightmm(:,3))];
    
    PS = polyshape(PSx,PSy);


end

function PS = analyseWorkspaceC(CentremmB,CentremmC, CentremmS2,CentremmS)

    Leftmm = CentremmS2;
    Rightmm = CentremmS;
    
    TipEnd_left = Leftmm(end,:);
    TipEnd_centre = CentremmB(end,:);
    TipEnd_right = Rightmm(end,:);
    
    
    % Combine the x and y coordinates into vectors
    x = [TipEnd_left(1), TipEnd_centre(1), TipEnd_right(1)];
    y = [TipEnd_left(3), TipEnd_centre(3), TipEnd_right(3)];
    
    % Fit a quadratic curve to the points
    p = polyfit(x, y, 2); 
    
    % Generate points along the fitted curve for plotting
    xFit = linspace(min(x), max(x), 101); % Generate 100 points for a smooth curve
    yFit = polyval(p, xFit); % Evaluate the polynomial at the xFit points
    
    Leftmm = CentremmC;
    
    PSx = [Leftmm(:,1);(xFit(51:end))';flip(Rightmm(:,1))];
    PSy = [Leftmm(:,3);(yFit(51:end))';flip(Rightmm(:,3))];
    
    PSx = [PSx;-PSx];
    PSy = [PSy;PSy];
    
    PS = polyshape(PSx,PSy);


end

%For dual material
function [PS,PSx,PSy] = tipWorkspace(CentremmB,CentremmC,CentremmC2,CentremmS,CentremmS2,CentremmJ,CentremmJ2)
    
    Leftmm = CentremmJ2;
    Rightmm = CentremmJ;
    
    TipEnd_left = Leftmm(end,:);
    TipEnd_centre = CentremmB(end,:);
    TipEnd_right = Rightmm(end,:);
    
    
    % Combine the x and y coordinates into vectors
    x = [TipEnd_left(1), TipEnd_centre(1), TipEnd_right(1)];
    y = [TipEnd_left(3), TipEnd_centre(3), TipEnd_right(3)];
    
    % Fit a quadratic curve to the points
    p = polyfit(x, y, 2); 
    
    % Get Curve for J to J
    xFit = linspace(min(x), max(x), 101); % Generate 100 points for a smooth curve
    yFit = polyval(p, xFit); % Evaluate the polynomial at the xFit points
    
    Leftmm = CentremmC;
     
    
    PSx = [flip(xFit');CentremmJ2(end,1)];
    PSy = [flip(yFit');CentremmJ2(end,3)];
    
    PS = polyshape(PSx,PSy);
end

% for double magnetic
function [PS,PSx,PSy] = tipWorkspace2(CentremmB,CentremmC,CentremmC2,CentremmS,CentremmS2,CentremmJ,CentremmJ2)
    
    Leftmm = CentremmS2;
    Rightmm = CentremmS;
    
    TipEnd_left = Leftmm(end,:);
    TipEnd_centre = CentremmB(end,:);
    TipEnd_right = Rightmm(end,:);
    
    
    % Combine the x and y coordinates into vectors
    x = [TipEnd_left(1), TipEnd_centre(1), TipEnd_right(1)];
    y = [TipEnd_left(3), TipEnd_centre(3), TipEnd_right(3)];
    
    % Fit a quadratic curve to the points
    p = polyfit(x, y, 2); 
    
    % Get Curve for J to J
    xFit = linspace(min(x), max(x), 101); % Generate 100 points for a smooth curve
    yFit = polyval(p, xFit); % Evaluate the polynomial at the xFit points
    
    %% Evaluate J curve intercept point
    curve2 = CentremmJ2(50:end,:);
    curve1 = CentremmJ(50:end,:);
    
    curve1 = sortrows(curve1, 1); % Sort by x values
    curve2 = sortrows(curve2, 1); % Sort by x values
    
    % Extract unique x-values from curve2
    [xUnique, ia, ic] = unique(curve2(:,1));
    
    % Average y-values for duplicate x-values in curve2
    yUnique = accumarray(ic, curve2(:,3), [], @mean);
    
    % Interpolate using unique x-values and their averaged y-values
    interpCurve2Y = interp1(xUnique, yUnique, curve1(:,1), 'linear', 'extrap');
    
    % Continue as before
    signChangeIndices = find(diff(sign(curve1(:,3) - interpCurve2Y)) ~= 0);
    xIntersect = (curve1(signChangeIndices, 1) + curve1(signChangeIndices+1, 1)) / 2;
    yIntersect = (curve1(signChangeIndices, 3) + interpCurve2Y(signChangeIndices+1)) / 2;

    % NOTE: specific to this picture only
    xIntersect = xIntersect(2);
    yIntersect = yIntersect(2);

    % Combine the x and y coordinates into vectors
    x2 = [CentremmC2(end,1), xIntersect, CentremmC(end,1)];
    y2 = [CentremmC2(end,3), yIntersect, CentremmC(end,3)];
    
    % Fit a quadratic curve to the points
    p2 = polyfit(x2, y2, 2); 
    
    % Get Curve for J to J
    xFit2 = linspace(min(x2), max(x2), 101); % Generate 100 points for a smooth curve
    yFit2 = polyval(p2, xFit2); % Evaluate the polynomial at the xFit points


    %% Evaluate PS
     
    
    PSx = [xFit';CentremmC(end,1);flip(xFit2');[CentremmC2(end,1),CentremmS2(end,1)]'];
    PSy = [yFit';CentremmC(end,3);flip(yFit2');[CentremmC2(end,3),CentremmS2(end,3)]'];
    
    PS = polyshape(PSx,PSy);

end

function img = combineImages(img1,img2,img3,img4)

    % Convert the image to grayscale if it is not already
    if size(img1, 3) == 3
        img1 = rgb2gray(img1);
    end

    % Convert the image to grayscale if it is not already
    if size(img2, 3) == 3
        img2 = rgb2gray(img2);
    end

    % Convert the image to grayscale if it is not already
    if size(img3, 3) == 3
        img3 = rgb2gray(img3);
    end

    % Convert the image to grayscale if it is not already
    if size(img4, 3) == 3
        img4 = rgb2gray(img4);
    end

    for i = 1:size(img1,1)
        for j = 1:size(img1,2)
            
            if (img1(i,j)==0 || img2(i,j) ==0 || img3(i,j)==0 || img4(i,j)==0)
                img(i,j) = 0;
            else
                img(i,j) = 255;
            end
        end
    end

end

function [X3D,Y3D,Z3D] = revolveTipWorkspace(PSx,PSz)

    % Convert the 2D profile to 3D by setting Y-coordinates to 0 (since it's in the X-Z plane)
    PSy = zeros(size(PSx)); % This step is more for clarity and understanding the transition to 3D
    
    % Define the angles for revolution
    angles = linspace(0, 2*pi, 360); % Full revolution around Z-axis


    
    % Initialize matrices to hold the coordinates of the 3D shape
    X3D = zeros(length(PSx), length(angles));
    Y3D = zeros(length(PSx), length(angles));
    Z3D = zeros(length(PSx), length(angles));

    
    % Calculate the coordinates for each point in the revolution
    for i = 1:length(PSx)
        X3D(i, :) = PSx(i) * cos(angles); % Calculate X-coordinates
        Y3D(i, :) = PSx(i) * sin(angles); % Calculate Y-coordinates (sine component gives the circular revolution)
        Z3D(i, :) = repmat(PSz(i), 1, length(angles)); % Z-coordinates remain constant for each point at all angles
    end
    

end