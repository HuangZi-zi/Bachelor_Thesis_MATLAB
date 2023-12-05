function [edgeImage] = u_find_edge(inputIMG)
%计算图像中的边缘


kernel1=strel('square',1);
kernel3=strel('square',3);
kernel5=strel('square',5);

% 获取原始图像的大小
[height, width, ~] = size(inputIMG);

% % Split the image into three color channels (Red, Green, Blue)
% redChannel = inputIMG(:, :, 1);
% greenChannel = inputIMG(:, :, 2);
blueChannel = inputIMG(:, :, 3);

% 使用 Sobel 边缘检测算法检测边缘
% edgeImage = edge(inputIMG, 'sobel');
% edgeRed = edge(redChannel, 'sobel');
% edgeGreen = edge(greenChannel, 'sobel');
edgeBlue = edge(blueChannel, 'sobel');

% Combine the three edge images into a single image
% combinedEdges = uint8(edgeRed+edgeGreen+edgeBlue) * 85;
% combinedEdges = imbinarize(combinedEdges,0.5);

% 裁剪边框
x = 10; % Starting column
y = 10;  % Starting row
w = width-2*x; % Width of the ROI
h = height-2*y; % Height of the ROI
% cropImage = combinedEdges(y:y+h-1, x:x+w-1);
cropImage = edgeBlue(y:y+h-1, x:x+w-1);

% 重新填充
paddedImage = zeros(height, width);
paddedImage(x + 1:x + h, y + 1:y + w) = cropImage;

edgeImage=paddedImage;

% % 开运算
% edgeImage = imerode(edgeImage,kernel1);
% edgeImage = imdilate(edgeImage,kernel1);
% 
% % 闭运算
% edgeImage = imdilate(edgeImage,kernel1);
% edgeImage = imerode(edgeImage,kernel1);


% % 模糊，使滑动窗口与边缘相交出现最值而非多个相同值
% edgeImage=double(edgeImage);
% edgeImageBlur=imgaussfilt(edgeImage, 5);
% 
% % 可视化结果
figure;
% subplot(2, 1, 1);
% imshow(edgeImage, 'InitialMagnification', 'fit');
% title('边缘检测结果');
% 
% subplot(2, 1, 2);
% imshow(edgeImageBlur);
% title('边缘检测输出');

% 可视化结果
% figure;
% subplot(2,2,1);
% imshow(cannyRed);
% title('Red');
% 
% subplot(2,2,2);
% imshow(cannyGreen);
% title('Green');
% 
% subplot(2,2,3);
% imshow(cannyBlue);
% title('Blue');

% subplot(2,2,4);
imshow(edgeImage);
title('Edges');
% figure;
% subplot(2, 1, 1);
% imshow(edgeImage, 'InitialMagnification', 'fit');
% title('边缘检测结果');