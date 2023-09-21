function [edgeImage] = u_find_edge(inputIMG)

kernel1=strel('square',1);
kernel3=strel('square',3);


% 获取原始图像的大小
[height, width] = size(inputIMG);

% 使用 Canny 边缘检测算法检测边缘
edgeImage = edge(inputIMG, 'Canny', 0.19);

% 裁剪边框
x = 10; % Starting column
y = 10;  % Starting row
w = width-2*x; % Width of the ROI
h = height-2*y; % Height of the ROI
cropImage = edgeImage(y:y+h-1, x:x+w-1, :);

% 重新填充
paddedImage = zeros(height, width);
paddedImage(x + 1:x + h, y + 1:y + w) = cropImage;


%腐蚀
edgeImage = imerode(paddedImage,kernel1);
% 膨胀，使边缘更加明显
edgeImage = imdilate(edgeImage,kernel3);

% % 模糊，使滑动窗口与边缘相交出现最值而非多个相同值
% edgeImage=double(edgeImage);
% edgeImageBlur=imgaussfilt(edgeImage, 5);
% 
% % 可视化结果
% figure;
% subplot(2, 1, 1);
% imshow(edgeImage, 'InitialMagnification', 'fit');
% title('边缘检测结果');
% 
% subplot(2, 1, 2);
% imshow(edgeImageBlur);
% title('边缘检测输出');

figure;
% subplot(2, 1, 1);
imshow(edgeImage, 'InitialMagnification', 'fit');
title('边缘检测结果');