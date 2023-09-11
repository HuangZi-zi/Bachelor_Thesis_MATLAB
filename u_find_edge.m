function [edgeImageBlur] = u_find_edge(inputIMG)

% 定义填充大小（根据需要调整）
paddingSize = 20; % 这里使用20像素的填充大小

% 获取原始图像的大小
[height, width] = size(inputIMG);

% 创建一个新的大图像，用于填充
paddedImage = zeros(height + 2 * paddingSize, width + 2 * paddingSize);
paddedImage(paddingSize + 1:paddingSize + height, paddingSize + 1:paddingSize + width) = inputIMG;

% 使用 Canny 边缘检测算法检测边缘
edgeImage = edge(paddedImage, 'Canny');

% 去除填充部分，保留原始图像部分
edgeImage = edgeImage(paddingSize + 1:paddingSize + height, paddingSize + 1:paddingSize + width);

% 膨胀，使边缘更加明显
kernel=strel('square',5);
edgeImage = imdilate(edgeImage,kernel);

% 模糊，使滑动窗口与边缘相交出现最值而非多个相同值
f_blur=fspecial("gaussian");
edgeImageBlur=imfilter(edgeImage,f_blur);

% 可视化结果
figure;
imshow(edgeImage, 'InitialMagnification', 'fit');
title('边缘检测结果');






