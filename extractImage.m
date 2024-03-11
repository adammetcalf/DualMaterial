close all;
clear;
clc;

%Load image
imgBase = imread('C:\Users\user\Desktop\Work\oppMagExperiment04032024\Base.png');
img = imread('C:\Users\user\Desktop\Work\oppMagExperiment04032024\2. FieldRtoL\Cshape.png');

% Rect = [xmin ymin width height]
rect = [1061.5-650,0,650*2,650*2];

imgBaseCropped = imcrop(imgBase,rect);
imgCropped = imcrop(img,rect);
clear('imgBase', 'img')

%figure(1)
%imshow(imgBaseCropped);

%figure(2)
%imshow(imgCropped);

%% Thresholding
grayImageBase = rgb2gray(imgBaseCropped);
grayImage = rgb2gray(imgCropped);

figure(3)
imshow(grayImageBase);

figure(4)
imshow(grayImage);

% Adaptive thresholding
bwAdaptiveBase = imbinarize(grayImageBase,'adaptive','ForegroundPolarity','dark','Sensitivity',0.4);
bwAdaptive = imbinarize(grayImage,'adaptive','ForegroundPolarity','dark','Sensitivity',0.4');

figure(5)
imshow(bwAdaptiveBase);

figure(6)
imshow(bwAdaptive);


%% Get image coordinates
%[x, y] = ginput(1);

%Display the coordinates
%disp(['You clicked at X: ', num2str(x), ', Y: ', num2str(y)]);