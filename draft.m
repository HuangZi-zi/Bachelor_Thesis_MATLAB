clear;
%img=imread("Resource\mwprinciple_webcam.png");
img=imread("Resource\corridor.jpg");
%tb_qrcodeDec_v2(img)

IMGgray=im2gray(img);
%模糊
IMGblur = imgaussfilt(IMGgray,2);
%膨胀
IMGdil=imdilate(IMGblur,strel('square',21));
%腐蚀
IMGero=imerode(IMGdil,strel('square',9));
% % 自适应二值化
% IMGbin=uint8(imbinarize(IMGero,'adaptive')*255);
IMGgray=uint8(IMGero);
% 中值滤波
IMGfil=medfilt2(IMGgray,[15,15]);
%腐蚀
IMGero=imerode(IMGfil,strel('square',5));
%膨胀
IMGdil=imdilate(IMGero,strel('square',5));

% Define the diameter of the circular mask (e.g., 37 pixels)
r = 18;                      %半径
[imgH,imgW] = size(IMGdil);
ac = linspace(0, 2*pi, 50);   %# approximate circle with 50 points




            




% Define your threshold 't'
t = 10;  % Adjust the threshold as needed

% Initialize a binary result image
result = zeros(imgH, imgW);
% result = zeros(size(image, 1), size(image, 2));
% Calculate the radius of the circular mask
% Iterate over the image
for i = 1:50:imgH
    for j = 1:50:imgW
        c=[j i];%圆心坐标
        %get circular mask
        mask = poly2mask(r*cos(ac)+c(1), r*sin(ac)+c(2), imgH, imgW);
        grayROI = immultiply(IMGdil,mask);
        imshow(grayROI);
        
        % Calculate the average gray level inside the circular mask
        average_roi = sum(grayROI(:))/(pi*r*r);
        
        % Calculate the gray level at the center of the mask circle
        center_gray_level = IMGdil(i, j);
        
        % Calculate the difference
        difference = abs(average_roi - center_gray_level);
        
        % Check if the difference is less than the threshold 't'
        if difference < t
            result(i, j) = 1;
        end
    end
end

imshow(result);


