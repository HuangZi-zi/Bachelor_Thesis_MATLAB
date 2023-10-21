% clear;
% %img=imread("Resource\mwprinciple_webcam.png");
% img=imread("Resource\corridor.jpg");
% %tb_qrcodeDec_v2(img)
%
% IMGgray=im2gray(img);
% %模糊
% IMGblur = imgaussfilt(IMGgray,2);
% %膨胀
% IMGdil=imdilate(IMGblur,strel('square',21));
% %腐蚀
% IMGero=imerode(IMGdil,strel('square',9));
% % % 自适应二值化
% % IMGbin=uint8(imbinarize(IMGero,'adaptive')*255);
% IMGgray=uint8(IMGero);
% % 中值滤波
% IMGfil=medfilt2(IMGgray,[15,15]);
% %腐蚀
% IMGero=imerode(IMGfil,strel('square',5));
% %膨胀
% IMGdil=imdilate(IMGero,strel('square',5));
%
%
%
%
%
% % Define the diameter of the circular mask (e.g., 37 pixels)
% r = 18;                      %半径
% [imgH,imgW] = size(IMGdil);
%
% % Define your threshold 't'
% t = 10;  % Adjust the threshold as needed
%
% % Initialize a binary result image
% result = zeros(imgH, imgW);
% % result = zeros(size(image, 1), size(image, 2));
% % Calculate the radius of the circular mask
% % Iterate over the image
% for i = 1:5:imgH
%     for j = 1:5:imgW
%         c=[j i];%圆心坐标
%         %get circular mask
%         % Extract the region of interest (ROI) using the circular mask
%         roi = IMGdil(max(1,(j-r)):min(imgH,(j+r)),max(1,(i-r)):min(imgW,(i+r)));
%
%         imshow(roi);
%
%         % Calculate the average gray level inside the circular mask
%         average_roi = sum(roi(:))/(4*r*r);
%
%         % Calculate the gray level at the center of the mask circle
%         center_gray_level = IMGdil(i, j);
%
%         % Calculate the difference
%         difference = abs(average_roi - center_gray_level);
%
%         % Check if the difference is less than the threshold 't'
%         if difference < t
%             result(i, j) = 1;
%         end
%     end
% end
%
% imshow(result);
%
%

clear obj2
obj2=serialport("COM9",19200,'Timeout', 0.2);
v=400;

%左转90度：vl=vr=200，t=5.9
%右转90度：vl=vr=200+32768, t=5.9

vl=v+32768;
vr=v+32768;
vlhex=dec2hex(vl,4);
vrhex=dec2hex(vr,4);
vlg= vlhex(1:2) ;%高位
vld=vlhex(3:4);%低位
vrg= vrhex(1:2) ;%高位
vrd=vrhex(3:4) ;%低位

sendbuff=zeros(1,9);
sendbuff(1)= hex2dec('55');
sendbuff(2)= hex2dec('aa');
sendbuff(3)= hex2dec('71');
sendbuff(4)= hex2dec('04');
sendbuff(5)= hex2dec('10');
sendbuff(6)= hex2dec(vlg);
sendbuff(7)= hex2dec(vld);
sendbuff(8)= hex2dec(vrg);
sendbuff(9)= hex2dec(vrd);
%校验和
%校验位=前面所有数据之和，取最后两位
add=sum(sendbuff,[1 2 3 4 5 6 7 8 9]);
two_bits_d=rem(add,256);%10进制下对256取余，在16进制下为2位
% two_bits_h= dec2hex(two_bits_d);% 发送数据以10进制存储，因此不需转换
sendbuff(10)= two_bits_d;
write(obj2,sendbuff,"uint8");

pause(2.95);

% 停止机器人
sendbuff(1)= hex2dec('55');
sendbuff(2)= hex2dec('aa');
sendbuff(3)= hex2dec('71');
sendbuff(4)= hex2dec('04');
sendbuff(5)= hex2dec('10');
sendbuff(6)= hex2dec('00');
sendbuff(7)= hex2dec('00');
sendbuff(8)= hex2dec('00');
sendbuff(9)= hex2dec('00');
sendbuff(10)= hex2dec('84');
write(obj2,sendbuff,"uint8")

%关闭串口
delete(obj2);
clear obj2;
