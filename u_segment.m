function [outputIMG] = u_segment(inputIMG)
originx=size(inputIMG,2);
originy=size(inputIMG,1);

% Specify the top-left corner (x, y) and the width and height of the region
x = 1; % Starting column
y = ceil(0.3*size(inputIMG,1));  % Starting row
width = size(inputIMG,2); % Width of the region
height = ceil(0.7*size(inputIMG,1)); % Height of the region

% Crop the image using indexing
outputIMG = inputIMG(y:min(y+height-1,originy), x:min(x+width-1,originx), :);

% figure;
% imshow(outputIMG);
% title('Cropped Image');

end