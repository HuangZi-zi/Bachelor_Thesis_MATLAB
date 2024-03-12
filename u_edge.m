function [img_edge] = u_edge(img)
%计算图像中的边缘

% 是否三个通道分别处理。0否，1是
div=1;

% 获取原始图像的大小
[h,w,channel]=size(img);

if channel == 3 && div == 1
    % Split the image into three color channels (Red, Green, Blue)
    redChannel = img(:, :, 1);
    greenChannel = img(:, :, 2);
    blueChannel = img(:, :, 3);
%     edgeRed = Sobel_filter(redChannel);
%     edgeGreen = Sobel_filter(greenChannel);
%     edgeBlue = Sobel_filter(blueChannel);
    edgeRed = edge(redChannel, 'sobel');
    edgeGreen = edge(greenChannel, 'sobel');
    edgeBlue = edge(blueChannel, 'sobel');
    % Combine the three edge images into a single image
    edgeImage = double(edgeRed+edgeGreen+edgeBlue);
    edgeImage = imbinarize(edgeImage,0.9);
    
%     figure;
%     subplot 131 
%     imshow(edgeRed);
%     subplot 132 
%     imshow(edgeGreen);
%     subplot 133 
%     imshow(edgeBlue);
elseif channel==3 && div==0
    gray=im2gray(img);
    edgeImage = edge(gray, 'sobel');
    edgeimage = Sobel_filter(gray);
    imshowpair(edgeImage,edgeimage,"montage");
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

% filted=medfilt2(paddedImage,[3,2]);
filted=ordfilt2(paddedImage,9,ones(3,3));
% filted=imerode(filted);
filted=ordfilt2(filted,9,ones(5,5));
img_edge=filted;

figure;
imshow(img_edge);
title('边缘');
end

%% 适于对角线检测的Sobel算子
% Sobel算子计算梯度得到边缘图像
function result_grad=Sobel_filter(img)
img=double(img);
[h,w,channel]=size(img);
% 先将原始灰度图向外扩展一圈,否则原始图像最外圈像素点得不到卷积结果
extended_image = zeros(h+2,w+2);
for i = 2:h+1
    for j = 2:w+1
        extended_image(i,j) = img(i-1,j-1);
    end
end
% 给出Sobel算子（对角线方向）
Sobel_x = [0, 1, 2;
          -1, 0, 1;
          -2, -1,0];
Sobel_y = Sobel_x';
% 初始化扩展图像的x,y方向梯度
gradx = zeros(h+2,w+2);
grady = zeros(h+2,w+2);
for i = 1:h
    for j = 1:w
        % current是原始图像当前被卷积的3×3矩阵
        current = extended_image(i:i+2,j:j+2);
        G_x = Sobel_x .* current;
        G_y = Sobel_y .* current;
        % 循环结束后gradx,grady即为扩展图像各像素点的x,y方向梯度值
        gradx(i+1,j+1) = sum(G_x(:));
        grady(i+1,j+1) = sum(G_y(:));
    end
end
% 将扩展一圈的图像复原(其实就是把扩展图像最外圈的0消去)
grad_x = zeros(h,w);
grad_y = zeros(h,w);
for i = 1:h
    for j = 1:w
        grad_x(i,j) = gradx(i+1,j+1);
        grad_y(i,j) = grady(i+1,j+1);
    end
end
% 求梯度相关参数值
grad = sqrt(grad_x.^2 + grad_y.^2);   % 得到图像的梯度幅度值
% 给函数输出的结果赋值
% result_grad = imbinarize(grad,0.1);
result_grad=uint8(grad);
end
