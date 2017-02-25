function [imout] = gauss2d(sigma, loc_x, loc_y, imsize)
% creates 2D gaussian with 1 std
    [rows, columns] = deal(imsize(1),imsize(2));
    backgroundGrayLevel = 1;
    windowSize = 150; % Could be random if you want.
    numberOfGaussians = 1;
    % Create one Gaussian.
    g = fspecial('gaussian', windowSize, sigma);
    grayImage = backgroundGrayLevel * ones(rows+2*windowSize, columns+2*windowSize);
    % Create random signs so that the Gaussians are
    % randomly brighter or darker than the background.
    s = 1;
    % Note: g and grayImage are floating point images, not uint8,
    % though you could modify the program to have them be uint8 if you wanted.
    % Get a list of random locations.
    centerRow = windowSize/2 + loc_x;
    centerCol = windowSize/2 + loc_y;
    % Place the Gaussians on the image at those random locations.
    grayImage(centerRow:centerRow+windowSize-1, centerCol:centerCol+windowSize-1) = ...
        grayImage(centerRow:centerRow+windowSize-1, centerCol:centerCol+windowSize-1) + ...
        s * g;
    
    imout = imcrop(grayImage,[windowSize windowSize columns-1 rows-1 ]);
end