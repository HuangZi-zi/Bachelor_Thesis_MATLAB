function [img_edge] = u_edge(img)
%计算图像中的边缘

% 是否三个通道分别处理。0否，1是
div=1;

% 获取原始图像的大小
[h,w,channel]=size(img);

if channel == 3 && div == 1
    % 将图像拆分成三个通道 (Red, Green, Blue)
    redChannel = img(:, :, 1);
    greenChannel = img(:, :, 2);
    blueChannel = img(:, :, 3);
    edgeRed = edge(redChannel, 'sobel');
    edgeGreen = edge(greenChannel, 'sobel');
    edgeBlue = edge(blueChannel, 'sobel');
    % 重新组合三个通道
    edgeImage = double(edgeRed+edgeGreen+edgeBlue);
    edgeImage = imbinarize(edgeImage,0.9);
elseif channel==3 && div==0
    gray=im2gray(img);
    edgeImage = edge(gray, 'sobel');
else
    edgeImage = edge(img, 'sobel');
end

% 裁剪边框
x = 10; % Starting column
y = 10;  % Starting row
w_roi = w-2*x; % Width of the ROI
h_roi = h-2*y; % Height of the ROI
cropImage = edgeImage(y:y+h_roi-1, x:x+w_roi-1);

% 重新填充
paddedImage = zeros(h, w);
paddedImage(x + 1:x + h_roi, y + 1:y + w_roi) = cropImage;

% 滤波
filted=ordfilt2(paddedImage,9,ones(3,3));
filted=ordfilt2(filted,9,ones(5,5));
img_edge=filted;

figure;
imshow(img_edge);
title('边缘');
end
