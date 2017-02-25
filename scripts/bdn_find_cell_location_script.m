% location registration

%% loc
% 250X150 pixels (example location [156 , 60]

[x,y] = deal(156,60);
%% registration frame

cd 'Z:\Noam\Code\RSgui\graphics'
bdn_template = 'Binary Dense Noise.png';
PosX = 315; % higher value -> lower location
PosY = 344; % higher value -> lefter location
heightScalers = 0.44;
[screenXpixels, screenYpixels] = Screen('WindowSize', window); % 800 x 600
[xCenter, yCenter] = RectCenter(windowRect);

% Get the aspect ratio of the image. We need this to maintain the aspect
% ratio of the image when we draw it different sizes. Otherwise, if we
% don't match the aspect ratio the image will appear warped / stretched
aspectRatio = s2 / s1;
imageHeights = screenYpixels .* heightScalers;
imageWidths = imageHeights .* aspectRatio;

% Make the destination rectangles for our image. We will draw the image
% multiple times over getting smaller on each iteration. So we need the big
% dstRects first followed by the progressively smaller ones
dstRects = zeros(4, numImages);
for i = 1:numImages
    theRect = [0 0 imageWidths(i) imageHeights(i)];
    dstRects(:, i) = CenterRectOnPointd(theRect, PosX,PosY);
end
Screen('DrawTextures', window, imageTexture, [], dstRects);
Screen('Flip', window);
%% frame 1
% see https://github.com/Psychtoolbox-3/Psychtoolbox-3/wiki/FAQ:-Screenshots-or-Recordings

imageArray = Screen('GetImage',window);
imwrite(imageArray, 'frame1.jpg')
%% stimulus mat frame (new bdn)
 
handles.stimulus_scale = '1';
set( handles.locationX,'String','400'); % doesn't matter
set( handles.locationY,'String','300');
cd 'C:\Users\noambox\Documents\Sync\stimulus files\src'
load('binary_dense_noise_stimulus_new.mat');
bdn_mat = M;
pixelSizes=Screen(whichScreen,'PixelSizes');
[window,screenRect] = Screen(whichScreen,'OpenWindow',0,[],max(pixelSizes));
[screenXpixels, screenYpixels] = Screen('WindowSize', window); % 800 x 600
aspectRatio = s2 / s1;
imageHeight = screenYpixels .* scale;
imageWidth = imageHeight .* aspectRatio;
% disp(['***** Intensity is: ', num2str(Intesity),' *******']);
theRect = [0 0 imageWidth imageHeight];
dstRects = CenterRectOnPointd(theRect, posX, posY); % image Position
stim = Screen('MakeTexture',window,bdn_mat(:,:,mod(k,stimMatSize)+1 ) );
Screen('DrawTexture',window,stim, [] ,dstRects ,[], [] ,[] ,[Intesity Intesity Intesity]);
%% stimulus mat frame (old bdn)

handles.stimulus_scale = '0.44';
%% frame 2

imageArray = Screen('GetImage',window);
imwrite(imageArray, 'frame2.jpg')
%% example with OneNote results 

x1 = 203; y1 = 56;
T1 = [[0.7156347255271192, -0.006083132214985391, 0]' [0.012243919072709428, -0.7209675926149547, 0]' [218.1354410707739, 402.5876157241306, 1]'];
T2 = [[0.5426369331949202, 6.417199507468235e-5, 0]' [5.326217994102055e-6, 0.5710082471585562, 0 ]' [72.71713363542452, 30.962022148317914, 1]'];
[w,h] = deal(200,150); % project to matrix dimensions
final = zeros(size(frame_mat));
[s,k] = deal(w/560, h/420);

coord_t = (T2 * (T1 * [x1, y1, 1]')); 
x_final = round( coord_t(1) * s ); 
y_final = round( coord_t(2) * k );
im = imread('C:\Users\noambox\Documents\Sync\code\Data Analysis\Visualization\bdn_location\frame1.jpg');
im_small = imresize(im , 1/4);
hf = figure; 
yellow = uint8([255 255 0]); % [R G B]; class of yellow must match class of I
shapeInserter = vision.ShapeInserter('Shape','Circles','BorderColor','Custom','CustomBorderColor',yellow);
circles = int32([x_final y_final 2]);
J = step(shapeInserter, im_small, circles);
imshow(J);
