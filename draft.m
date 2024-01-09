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

% clear obj2
% obj2=serialport("COM9",19200,'Timeout', 0.2);
% v=400;
% 
% %左转90度：vl=vr=200，t=5.9
% %右转90度：vl=vr=200+32768, t=5.9
% 
% vl=v+32768;
% vr=v+32768;
% vlhex=dec2hex(vl,4);
% vrhex=dec2hex(vr,4);
% vlg= vlhex(1:2) ;%高位
% vld=vlhex(3:4);%低位
% vrg= vrhex(1:2) ;%高位
% vrd=vrhex(3:4) ;%低位
% 
% sendbuff=zeros(1,9);
% sendbuff(1)= hex2dec('55');
% sendbuff(2)= hex2dec('aa');
% sendbuff(3)= hex2dec('71');
% sendbuff(4)= hex2dec('04');
% sendbuff(5)= hex2dec('10');
% sendbuff(6)= hex2dec(vlg);
% sendbuff(7)= hex2dec(vld);
% sendbuff(8)= hex2dec(vrg);
% sendbuff(9)= hex2dec(vrd);
% %校验和
% %校验位=前面所有数据之和，取最后两位
% add=sum(sendbuff,[1 2 3 4 5 6 7 8 9]);
% two_bits_d=rem(add,256);%10进制下对256取余，在16进制下为2位
% % two_bits_h= dec2hex(two_bits_d);% 发送数据以10进制存储，因此不需转换
% sendbuff(10)= two_bits_d;
% write(obj2,sendbuff,"uint8");
% 
% pause(2.95);
% 
% % 停止机器人
% sendbuff(1)= hex2dec('55');
% sendbuff(2)= hex2dec('aa');
% sendbuff(3)= hex2dec('71');
% sendbuff(4)= hex2dec('04');
% sendbuff(5)= hex2dec('10');
% sendbuff(6)= hex2dec('00');
% sendbuff(7)= hex2dec('00');
% sendbuff(8)= hex2dec('00');
% sendbuff(9)= hex2dec('00');
% sendbuff(10)= hex2dec('84');
% write(obj2,sendbuff,"uint8")
% 
% %关闭串口
% delete(obj2);
% clear obj2;

%% 直线检测
% img=imread("Resource\conor.png");
% % img =imcrop(img,[90,450,1100,400]);
% % img = u_basic_process(img);
% edge_img=u_edge(img);
% laneLines=u_line_hough(img, edge_img);
% points=u_point_orb(img);
% %x=ay+b
% % y=1:size(canny,1);
% % x_mid=y.*(laneLines(2).a+laneLines(1).a)/2+(laneLines(2).b+laneLines(1).b)/2;
% % hold on
% % plot(x_mid,y,'.');
% 
% u_ROI(laneLines,img);

%% 角点检测
% close all;
% clc;
% % 读取图像信息（原图为灰度图）
% img = imread("Resource\conor.png");
% blueChannel = img(:, :, 3);
% img =im2gray(blueChannel);
% img = imgaussfilt(img,2);
% 
% [m,n] = size(img);
% % 先在原图外围扩展一圈
% tmp = zeros(m+2,n+2);
% tmp(2:m+1,2:n+1) = img;
% % 初始化各一阶偏导矩阵
% Ix = zeros(m+2,n+2);
% Iy = zeros(m+2,n+2);
% E = zeros(m+2,n+2);
% % 求偏导
% Ix(:,2:n) = tmp(:,3:n+1) - tmp(:,1:n-1);
% Iy(2:m,:) = tmp(3:m+1,:) - tmp(1:m-1,:);
% Ix2 = Ix(2:m+1,2:n+1) .^ 2;
% Iy2 = Iy(2:m+1,2:n+1) .^ 2;
% Ixy = Ix(2:m+1,2:n+1) .* Iy(2:m+1,2:n+1);
% %生成高斯卷积核，便于对Ix2、Iy2、Ixy进行平滑
% % sigma = 2
% h = fspecial('gaussian',[3 3],2);
% Ix2 = filter2(h,Ix2);
% Iy2 = filter2(h,Iy2);
% Ixy = filter2(h,Ixy);
% 
% % 初始化Rmax
% Rmax = 0;
% R = zeros(m,n);
% for i = 1 : m
%     for j = 1 : n
%         M = [Ix2(i,j) Ixy(i,j);
%              Ixy(i,j) Iy2(i,j)];
%         R(i,j) = det(M) - 0.06 * (trace(M))^2;
%         if R(i,j) > Rmax
%             Rmax = R(i,j);
%         end
%     end
% end
% % 显示图像
% imshow(img);
% title('角点检测');
% hold on;
% 
% % 求角点
% tmp(2:m+1,2:n+1) = R;
% result = zeros(m+2,n+2);
% result(2:m+1,2:n+1) = img;
% for i = 2 : m + 1
%     for j = 2 : n + 1
%         % 阈值为0.02*Rmax
%         % 认为R值大于阈值的点为角点
%         % 求当前像素点的邻域
%         current = [tmp(i-1,j-1) tmp(i-1,j) tmp(i-1,j+1);
%                    tmp(i,j-1)   tmp(i,j)   tmp(i,j+1);
%                    tmp(i+1,j-1) tmp(i+1,j) tmp(i+1,j+1)];
%         % 若当前像素点的R值大于阈值且它是其八邻域中R值最大的点，则它为角点
%         if tmp(i,j) >= 0.02 * Rmax && tmp(i,j) >= max(max(current))
%                 result(i,j) = 255;
%                 % plot绘制点的时候是以左上角为原点，水平向右为x正半轴轴，竖直向下为y正半轴
%                 % 这和我们对于图像矩阵坐标的直观印象恰好相反
%                 plot(j,i,'go')
%         end
%     end
% end
% 
% % 这是测试plot绘制点的代码
% % for i = 1 : m
% %     for j = 1 : n
% %         plot(i,j,'b+');
% %         pause;
% %     end
% % end

%% 
% 
% rng('default') % For reproducibility
% X = rand(3,2);
% D2 = pdist(X,@naneucdist);
% 
% function D2 = naneucdist(XI,XJ)  
% %NANEUCDIST Euclidean distance ignoring coordinates with NaNs
% n = size(XI,2);
% sqdx = (XI-XJ).^2;
% nstar = sum(~isnan(sqdx),2); % Number of pairs that do not contain NaNs
% nstar(nstar == 0) = NaN; % To return NaN if all pairs include NaNs
% D2squared = sum(sqdx,2,'omitnan').*n./nstar; % Correction for missing coordinates
% D2 = sqrt(D2squared);
% end

%% 
% rng('default') % For reproducibility
% X = rand(10,3).*10;
% tree = linkage(X,'average');
% c=cluster(tree,'MaxClust',5);
% figure()
% scatter3(X(:,1),X(:,2),X(:,3),10,c)
% figure()
% dendrogram(tree)

%% 调用kinect

addpath('C:\Users\YawnFun\Documents\project\Kin2\Mex');
clear
close all

% Create Kinect 2 object and initialize it
% Available sources: 'color', 'depth', 'infrared', 'body_index', 'body',
% 'face' and 'HDface'
k2 = Kin2('color','depth');
i=1;
% images sizes
depth_width = 512; depth_height = 424; outOfRange = 4096;
color_width = 1920; color_height = 1080;

% Color image is to big, let's scale it down
colorScale = 1;

% Create matrices for the images
depth = zeros(depth_height,depth_width,'uint16');
color = zeros(color_height*colorScale,color_width*colorScale,3,'uint8');
depthColor = zeros(depth_height,depth_width,3,'uint8');

% depth stream figure
% figure, h1 = imshow(depth,[0 outOfRange]);
% title('Depth Source (press q to exit)')
% colormap('Jet')
% colorbar
% set(gcf,'keypress','k=get(gcf,''currentchar'');'); % listen keypress

% color stream figure
figure, h2 = imshow(color,[]);
title('Color Source (press q to exit)');
set(gcf,'keypress','k=get(gcf,''currentchar'');'); % listen keypress

% depthColor stream figure
figure, h4 = imshow(depthColor,[]);
title('Color2DephAligned (press q to exit)');
set(gcf,'keypress','k=get(gcf,''currentchar'');'); % listen keypress

% Loop until pressing 'q' on any figure
k=[];

disp('Press q on any figure to exit')

while true
    % Get frames from Kinect and save them on underlying buffer
    validData = k2.updateData;
    
    % Before processing the data, we need to make sure that a valid
    % frame was acquired.
    if validData
        % Copy data to Matlab matrices        
        depth = k2.getDepth;
        color = k2.getColor;
        depthColor = k2.getAlignColor2Depth;               
        
        % update depth figure
%         depth(depth>outOfRange) = outOfRange; % truncate depht
%         set(h1,'CData',depth); 

        % update color figure
        color = imresize(color,colorScale);
        set(h2,'CData',color); 
       
        set(h4,'CData',depthColor);
        % u_plane_regiongrowing(color,depth);
    end
    
    
    % If user presses 'q', exit loop
    if ~isempty(k)
        if strcmp(k,'s')
            imwrite(color,"Resource\snapc"+i+".png");
            imwrite(depthColor,"Resource\snapcd"+i+".png");
            imwrite(depth.*16,"Resource\snapd"+i+".png");

            disp("snapshot");
            i=i+1;
            k=[];
        elseif strcmp(k,'q')
            break
        end
        
    end
  
    pause(0.02)
end

% Close kinect object
k2.delete;

close all;

%% 测试坐标系的相互映射关系
addpath('C:\Users\YawnFun\Documents\project\Kin2\Mex');
clear
close all
k2 = Kin2('color','depth');

x=(1:600)';
y=(1:600)';
% x1=[1;512];
% y1=[1;424];
% y=ones(30,1).*20;

% ans1=k2.mapDepthPoints2Color([x,y]);
ans2=k2.mapColorPoints2Depth([300,800]);