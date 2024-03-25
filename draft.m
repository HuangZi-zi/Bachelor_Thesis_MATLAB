%% 串口收发
    clear obj2
    obj2=serialport("COM3",19200,'Timeout', 0.2);
    sendbuff=zeros(1,9);
    sendbuff(1)= hex2dec('55');
    sendbuff(2)= hex2dec('aa');
    sendbuff(3)= hex2dec('71');
    sendbuff(4)= hex2dec('04');
    sendbuff(5)= hex2dec('10');
    sendbuff(6)= hex2dec('03');
    sendbuff(7)= hex2dec('20');
    sendbuff(8)= hex2dec('03');
    sendbuff(9)= hex2dec('20');
    sendbuff(10)= hex2dec('ca');
    write(obj2,sendbuff,"uint8")

%% 论证使用RGB处理的必要性

img=imread("Resource\mainline.jpg");
imggray=im2gray(img);
% imshowpair(img,imggray,"montage");

u_edge(img);
%%
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
img=imread("Resource\corridor.jpg");
% img =imcrop(img,[90,450,1100,400]);
% img = u_basic_process(img);
edge_img=u_edge(img);
laneLines=u_line_hough(img, edge_img);
% points=u_point_orb(img);
%x=ay+b
% y=1:size(canny,1);
% x_mid=y.*(laneLines(2).a+laneLines(1).a)/2+(laneLines(2).b+laneLines(1).b)/2;
% hold on
% plot(x_mid,y,'.');

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

% Create Kinect 2 object and initiqalize it
% Available sources: 'color', 'depth', 'infrared', 'body_index', 'body',
% 'face' and 'HDface'
k2 = Kin2('color','depth');
i=30;
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

%% 测试点云处理平面提取
addpath('C:\Users\YawnFun\Documents\project\Kin2\Mex');
clear
close all

k2 = Kin2('color','depth');

% images sizes
depth_width = 512; depth_height = 424; outOfRange = 4000;

% Create matrices for the images
depth = zeros(depth_height,depth_width,'uint16');
pc = pointCloud(zeros(depth_height*depth_width,3));

% depth stream figure
% figure, h1 = imshow(depth,[0 outOfRange]);
% title('Depth Source (close figure to exit)')
% colormap('Jet')
% colorbar

% point cloud figure
pcFig.h = figure;
pcFig.ax = pcshow(pc);

disp('Close any figure to exit')
downsample = 2; % subsample pointcloud

% Main Loop
while true
    tic;
    % Get frames from Kinect and save them on underlying buffer
    validData = k2.updateData;
    
    % Before processing the data, we need to make sure that a valid
    % frame was acquired.
    if validData
        % Copy data to Matlab matrices
        depth = k2.getDepth;

        % update depth figure
        depth(depth>outOfRange) = outOfRange; % truncate depht
        
        % Display the depth image, 
        % if the user closes the window, the program ends
%         try
%             set(h1,'CData',depth); 
%         catch
%             break; % break the main loop 
%         end
          
        % Get the pointcloud with color from the Kinect
        % Select the output 'pointCloud' to use the MATLAB built-in
        % pointCloud object. 
        % For MATLAB versions older than 2015b, use 'output','raw' and use
        % scatter3 to plot the point cloud. See pointCloudDemo1.m
        pc = k2.getPointCloud('output','pointCloud','color','true');
        
        % Display the point cloud,
        % if the user closes the window, the program ends
        try
%             pcshow(pc,'Parent',pcFig.ax,'VerticalAxis','Y');
%             title(pcFig.ax,'Point Cloud');
%             xlabel(pcFig.ax,'X'); ylabel(pcFig.ax,'Y'); zlabel(pcFig.ax,'Z');
%             axis(pcFig.ax,[-4 4 -4 4 -4 4]);
            minDistance = 0.5;
            [labels,numClusters] = pcsegdist(pc,minDistance);
            pcshow(pc.Location,labels)
            colormap(hsv(numClusters))
            title('Point Cloud Clusters')
        catch
            break; % break the main loop
        end
    end
    toc;
    pause(0.02);
end
% Close kinect object
k2.delete;

%% 使用卷积的方法得到角点
img=imread("Resource\conor.png");
kernal=[0,0,0;1,1,1;-1,1,-1]';
gray=im2gray(img);
out=conv2(kernal,gray);
out2=imbinarize(out,300);
imshowpair(out2,gray,"montage");


%% SIFT手工编程模板
%使用SIFT为模板修改而来
I=imread("Resource\conor.png");
I=rgb2gray(I);

%sift算法提取特征位置，并给出描述子

pic_orig=I;
figure(1);
imshow(pic_orig);%imshow函数用来显示原始标准影像
title('原始标准影像');
 
%判断是否为彩色图像，是彩色则转换为灰度图像
[x,y,z]=size(pic_orig);
if(z~=1)
    pic_gray = rgb2gray(pic_orig);%rgb2gray将彩色图像转化为灰度图像
else
    pic_gray = pic_orig;
end
figure(2);
imshow(pic_gray);
title('灰度图像');


%构建高斯金字塔
[M,N]=size(pic_gray);
O=floor(log2(min(M,N)))-3;%确定高斯金字塔的组数(floor函数用来取整，round函数用来四舍五入)
S=3;%取高斯差分金字塔每组的层数为3
sigama0=1.6;%初始sigama0设为1.6
k=2^(1/S);
pic_gray0=imresize(pic_gray,2,'nearest');%先将原图像长、宽各扩展一倍
for i=1:O+1%遍历每一组
    if i==1
        pic_gray0=pic_gray0;
    else
        pic_gray0=imresize(gauss_pyr_img{i-1}{S+1},0.5,'nearest');%下一组第一层等于上一组倒数第三层直接降采样得到
    end
    for j=1:S+3%对每一组的每一层卷积
        sigama1=sigama0*2^(i-1);
        sigama=sigama1*k^(j-1);
        gauss_pyr_img{i}{j}=gauss(pic_gray0,sigama);
    end
end
 
%构建高斯差分金字塔
for i=1:O+1
    for j=1:S+2
        FOGDD{i}{j}=gauss_pyr_img{i}{j+1}-gauss_pyr_img{i}{j};
    end
end

for i=1:O+1
    for j=1:S+1
        SOGDD{i}{j}=FOGDD{i}{j+1}-FOGDD{i}{j};
    end
end


%DoG空间极值探测并去除低对比度的边缘相应点
location=cell(size(DoG_pyr_img));%用来存储极值点
curvature_threshold=10.0;
curvature_threshold=(curvature_threshold+1)^2/curvature_threshold;%主曲率阈值
contrast_threshold=0.03;%对比度阈值
xx=[1 -2 1];
yy=xx';
xy=[-1 0 1;0 0 0;1 0 -1]/4;
for i=1:O+1
    [M,N]=size(DoG_pyr_img{i}{1});
    for j=2:S+1
        contrast_mask=abs(DoG_pyr_img{i}{j}(:,:))>=contrast_threshold;%对比度高于阈值的像素值为1，低于的为0
        location{i}{j-1}=zeros(M,N);
        for ii=2:M-1
            for jj=2:N-1
                if(contrast_mask(ii,jj)==1)%筛选出对比度高于阈值的点
                    tmp1=DoG_pyr_img{i}{j-1}((ii-1):(ii+1),(jj-1):(jj+1)); 
                    tmp2=DoG_pyr_img{i}{j}((ii-1):(ii+1),(jj-1):(jj+1));
                    tmp3=DoG_pyr_img{i}{j+1}((ii-1):(ii+1),(jj-1):(jj+1));
                    tmp=[tmp1;tmp2;tmp3];
                    center=tmp2(2,2);%中心点
                    if(center==max(tmp(:))||center==min(tmp(:)))
                        Dxx=sum(tmp2(2,1:3).*xx);
                        Dyy=sum(tmp2(1:3,2).*yy);
                        Dxy=sum(sum(tmp2(:,:).*xy));%计算极值处的Hessian矩阵
                        trH=Dxx+Dyy;
                        detH=Dxx*Dyy-(Dxy)^2;
                        curvature=(trH^2)/detH;%计算主曲率
                        if(curvature<curvature_threshold)%主曲率小于阈值
                            location{i}{j-1}(ii,jj)=255;
                        end
                    end
                end
            end
        end
    end
end
 
%求取特征点的主方向
length=cell(size(location));%用来存储关键点梯度大小
direction=cell(size(location));%用来存储关键点梯度方向
sigama0=1.6;
count=0;%设置一个计数器，用来确定计算以特征点为中心图像幅角和幅值的区域半径
hist=zeros(1,36);%存储邻域内梯度方向直方图
for i=1:O+1
    [M,N]=size(gauss_pyr_img{i}{1});
    for j=2:S+1
        count=count+1;
        sigama=sigama0*k^count;%确定区域半径的sigama
        length{i}{j-1}=zeros(M,N);
        direction{i}{j-1}=zeros(M,N);
        r=8;%区域半径设为8
        for ii=r+2:M-r-1
            for jj=r+2:N-r-1
                if(location{i}{j-1}(ii,jj)==255)
                    for iii=ii-r:ii+r
                        for jjj=jj-r:jj+r
                            m=sqrt((gauss_pyr_img{i}{j}(iii+1,jjj)-gauss_pyr_img{i}{j}(iii-1,jjj))^2+(gauss_pyr_img{i}{j}(iii,jjj+1)-gauss_pyr_img{i}{j}(iii,jjj-1))^2);
                            theta=atan((gauss_pyr_img{i}{j}(iii,jjj+1)-gauss_pyr_img{i}{j}(iii,jjj-1))/(gauss_pyr_img{i}{j}(iii+1,jjj)-gauss_pyr_img{i}{j}(iii-1,jjj)));
                            theta=theta/pi*180;%将弧度化为角度
                            if(theta<0)
                                theta=theta+360;
                            end
                            w=exp(-(iii^2+jjj^2)/(2*(1.5*sigama)^2));%生成邻域各像元的高斯权重
                            if(isnan(theta))
                                if(gauss_pyr_img{i}{j}(iii,jjj+1)-gauss_pyr_img{i}{j}(iii,jjj-1)>=0)
                                    theta=90;
                                else
                                    theta=270;
                                end
                            end
                            t=floor(theta/10)+1;
                            hist(t)=hist(t)+m*w;%将幅值*高斯权重加入对应的直方图中
                        end
                    end
                    for iiii=2:35
                        hist(iiii)=(hist(iiii-1)+hist(iiii)+hist(iiii+1))/3;%直方图的平滑处理
                    end
                    for iiii=2:35
                        hist(iiii)=(hist(iiii-1)+hist(iiii)+hist(iiii+1))/3;%第二次平滑处理
                    end
                    [u,v]=max(hist(:));
                    length{i}{j-1}(ii,jj)=u;%存储主方向上的幅值*高斯权重
                    direction{i}{j-1}(ii,jj)=(v-1)*10;%存储主方向
                end
            end
        end
    end
end
                    
% %特征描述符的生成
% sigama0=1.6;
% count=0;
% description=[];%用来存储描述子
% index=0;%用来索引描述子
% description_1=cell(size(location));%用来存储索引，通过索引找到description中对应的描述子
% d=[];%存储邻域内每个像素梯度方向
% l=[];%存储邻域内每个像素梯度幅值
% f=zeros(1,8);%用来存放4*4邻域内8维描述子
% description_2=[];%用来存放128维描述子
% aaa=[];
% for i=1:O+1
%     [M,N]=size(gauss_pyr_img{i}{1});
%     for j=2:S+1
%         description_1{i}{j-1}=zeros(M,N);
%         count=count+1;
%         sigama=sigama0*k^count;
%         %r=floor((3*sigama*sqrt(2)*5+1)/2);%确定描述子所需要的图像区域半径
%         r=8;%设描述子所需要的图像区域半径
%         for ii=r+2:M-r-1
%             for jj=r+2:N-r-1
%                 if(length{i}{j-1}(ii,jj)~=0)
%                     theta_1=direction{i}{j-1}(ii,jj);%该邻域的主方向
%                     index=index+1;
%                     description_2=[];
%                     d=[];
%                     l=[];
%                     for iii=[ii-r:1:ii-1,ii+1:1:ii+r]
%                         for jjj=[jj-r:1:jj-1,jj+1:1:jj+r]
%                             m=sqrt((gauss_pyr_img{i}{j}(iii+1,jjj)-gauss_pyr_img{i}{j}(iii-1,jjj))^2+(gauss_pyr_img{i}{j}(iii,jjj+1)-gauss_pyr_img{i}{j}(iii,jjj-1))^2);
%                             theta=atan((gauss_pyr_img{i}{j}(iii,jjj+1)-gauss_pyr_img{i}{j}(iii,jjj-1))/(gauss_pyr_img{i}{j}(iii+1,jjj)-gauss_pyr_img{i}{j}(iii-1,jjj)));
%                             theta=theta/pi*180;%将弧度化为角度
%                             if(theta<0)
%                                 theta=theta+360;
%                             end
%                             w=exp(-(iii^2+jjj^2)/(2*(1.5*sigama)^2));%生成邻域各像元的高斯权重
%                             if(isnan(theta))
%                                 if(gauss_pyr_img{i}{j}(iii,jjj+1)-gauss_pyr_img{i}{j}(iii,jjj-1)>=0)
%                                     theta=90;
%                                 else
%                                     theta=270;
%                                 end
%                             end
%                             theta=theta+360-theta_1;%逆时针旋转至主方向
%                             theta=mod(theta,360);%取余
%                             d=[d,theta];
%                             l=[l,m*w];
%                         end
%                     end
%                     d=reshape(d,16,16);
%                     l=reshape(l,16,16);%将一维数组变为二维矩阵
%                     for r1=1:4
%                         for c1=1:4
%                             for iiii=1+(r1-1)*4:4*r1
%                                 for jjjj=1+(c1-1)*4:4*c1
%                                     t=floor(d(iiii,jjjj)/45+1);%方向
%                                     f(t)=f(t)+l(iiii,jjjj);
%                                 end
%                             end
%                             description_2=[description_2,f(:)];%得到一个128维的描述子
%                         end
%                     end
%                     description_2=description_2./sqrt(sum(description_2(:)));%归一化处理
%                     description=[description;description_2(:)];
%                     description_1{i}{j-1}(ii,jj)=index;
%                     aaa=[aaa;ii,jj];
%                 end
%             end
%         end
%     end
% end
% description=reshape(description,[],128);%将描述子变为128维
 
%绘制特征点
figure(3)
imshow(pic_gray)
hold on
plot(aaa(:,1),aaa(:,2),'y+')
%description_1{i}{j}(ii,jj)中不为0的点即为特征点，其值为一个索引号，用这个索引号可以到description中找到相应的描述子

% function y=gauss(x,sigama)
% %y--滤波后的图像
% %x--需滤波的图像
% %sigama--高斯滤波核；
% n=round(6*sigama);%四舍五入
% if mod(n,2)==0
%     n=n+1;%确定高斯卷积核的大小
% end
% h=fspecial('gaussian',[n,n],sigama);%创建卷积窗口
% y=filter2(h,x);%filter2函数进行二维图像的卷积
% end

%% FOGDD&SOGDD角点判别
% 旋转矩阵的计算
theta=0:45:360;
R_k=[cos(theta),sin(theta);-sin(theta),cos(theta)];

% function  [g,phi,psi]=gauss(n,sigma)
% % g 高斯滤波器；n坐标[x, y]'；sigma 方差
% g=e^(n'*n/(2*sigma^2))/(2*pi*sigma^2);
% phi=-([cos(theta),sin(theta)]*n*g)/sigma^2;
% end

%% 将深度与彩色对齐图像的基本处理
imgcd=imread("Resource\snapcd25.png");
imgd=imread("Resource\snapd25.png");

R=imgcd(:,:,1);
G=imgcd(:,:,2);
B=imgcd(:,:,3);
figure();
imshow(R);

% 利用最小值滤波得到高光遮罩，并去除高光
min_R = ordfilt2(R,1,ones(5,5));
min_G = ordfilt2(G,1,ones(5,5));
min_B = ordfilt2(B,1,ones(5,5));

ave_R=mean(min_R,"all");
ave_G=mean(min_G,"all");
ave_B=mean(min_B,"all");

min_R(min_R(:)<ave_R)=0;
min_G(min_G(:)<ave_G)=0;
min_B(min_B(:)<ave_B)=0;

figure();imshow(min_R);

R1=(R-0.3.*min_R);
G1=(R-0.3.*min_G);
B1=(R-0.3.*min_B);

% 孔洞填充
R1=imopen(R1,ones(3,3));
R1=imfill(R1);
% R1=imclose(R1,ones(3,3));
% figure();imshow(R1);
G1=imopen(G1,ones(3,3));
G1=imfill(G1);
B1=imopen(B1,ones(3,3));
B1=imfill(B1);

% 中值滤波
R1=ordfilt2(R1,1,ones(3,3));% 输入图像，序号，邻域。程序会对邻域内像素升序排列，将序号对应的值赋予中心像素
G1=ordfilt2(G1,1,ones(3,3));
B1=ordfilt2(B1,1,ones(3,3));

img(:,:,1)=R1;
img(:,:,2)=G1;
img(:,:,3)=B1;
figure();imshow(img);

u_plane_regiongrowing(imgcd,imgd);
u_plane_regiongrowing(img,imgd)
% 

%% 对比一次二维中值和两次一维中值
img=imread("Resource\lena.png");
img=im2gray(img);
subplot 221;imshow(img);title("origin");
img=imnoise(img,'gaussian',0.002);
img=imnoise(img,'salt & pepper',0.02);
img=im2double(img);
h=480;
w=640;

subplot 222;imshow(img);title("noise");
% 中值滤波
filter_size=5;
% hist_eq_double = im2double(hist_eq);
% x_dir=hist_eq_double(:)';
% tic;
x_dir=img(:)';
fil_x=medfilt1(x_dir,filter_size);
fil_x_re=reshape(fil_x,h,w)';
y_dir=fil_x_re(:)';
fil_xy=medfilt1(y_dir,filter_size);%2次1维中值0.03024628s
fil_xy_re=reshape(fil_xy,w,h)';
% toc;
% disp(['运行时间: ',num2str(toc)]);
tic
filted=medfilt2(img,[filter_size,filter_size]);%1次2维中值0.18264729s
toc
disp(['运行时间: ',num2str(toc)]);
subplot 223;imshow(fil_xy_re);title("xy one-dimension filter");
subplot 224;imshow(filted);title("two-dimensions filter");

%% 对比sobel和canny算子
img=imread("Resource\lena.png");
img=im2gray(img);
% subplot 221;imshow(img);title("origin");
img=imnoise(img,'gaussian',0.002);
img=imnoise(img,'salt & pepper',0.02);
img=im2double(img);
h=480;
w=640;

subplot 221;imshow(img);title("noise");
% 中值滤波
filter_size=5;
% hist_eq_double = im2double(hist_eq);
% x_dir=hist_eq_double(:)';
% tic;
x_dir=img(:)';
fil_x=medfilt1(x_dir,filter_size);
fil_x_re=reshape(fil_x,h,w)';
y_dir=fil_x_re(:)';
fil_xy=medfilt1(y_dir,filter_size);%2次1维中值0.03024628s
fil_xy_re=reshape(fil_xy,w,h)';
subplot 222;imshow(fil_xy_re);title("filtered");

edgesobel = edge(fil_xy_re, 'sobel');
edgecanny = edge(fil_xy_re, "canny");
subplot 223;imshow(edgesobel);title("edge-Sobel");
subplot 224;imshow(edgecanny);title("edge-Canny");

%% 曲线测试情况
img=imread("Resource\snapc25.png");
% edge_img=u_edge(img);
% u_line_hough(img,edge_img);
u_plane_regiongrowing(img,img);

%% 人工势场法
img=imread("Resource\snapc30.png");
edges=u_plane_regiongrowing(img,img);
[out,dir]=u_APF(img,edges);
imshow(out);
disp(dir);

%% 深度与彩色对齐
color=imread("Resource\snapc30.png");
depth=imread("Resource\snapd30.png");
color=imresize(color,[375,667]);
depthColor_c=fliplr(imcrop(color,[89 1 511 375]));
depthColor_d=fliplr(imcrop(depth,[1 8 511 374]));
% depthColor_c=im2gray(depthColor_c);
% imshowpair(depthColor_d,depthColor_c,"blend");
edges=u_plane_regiongrowing(depthColor_c,depthColor_d);
[out,dir]=u_APF(depthColor_c,edges);
figure(1);imshow(out);
disp(dir);

%% 新的邻域生长法编程
color=imread("Resource\snapc30.png");
depth=imread("Resource\snapd30.png");
color=imresize(color,[375,667]);
img_color=fliplr(imcrop(color,[89 1 511 375]));
img_depth=fliplr(imcrop(depth,[1 8 511 374]));

[height,width,channel]=size(img_color);
% 720,1080
nodesize=9;
if rem(height,nodesize)% 如果不能整除，则进一步裁剪
    img_color=imcrop(img_color,[1 1 width-1 nodesize*floor(height/nodesize)-1]);
end
if rem(width,nodesize)
    img_color=imcrop(img_color,[1 1 nodesize*floor(width/nodesize)-1 height-1]);
end
[height,width,channel]=size(img_color);

node_height=height/nodesize;
node_width=width/nodesize;
node_mid=floor(node_width/2);

nodes=cell(node_height, node_width);% nodes
depth=zeros(node_height, node_width);
edges=ones(node_height,2).*node_mid;
scan_lines=cell(node_height, node_width);
% node_last_row=cell(1, node_n);

% 划分为node
if rem(nodesize,2) % nodesize为奇数，取中间1行
    for i=1:node_height
        for j=1:node_width
            nodes{i,j}=img_color((i-1)*nodesize+1:min(i*nodesize,height), ...
                                 (j-1)*nodesize+1:min(j*nodesize,width),:);
            depth(i,j)=mean(img_depth((i-1)*nodesize+1:min(i*nodesize,height), ...
                                      (j-1)*nodesize+1:min(j*nodesize,width)),"all");
            scan_lines{i,j}=[(j-1)*nodesize+1:j*nodesize;...
                            repmat((i-1)*nodesize+1+floor(nodesize/2),1,nodesize);...
                            img_depth((i-1)*nodesize+1+floor(nodesize/2),(j-1)*nodesize+1:j*nodesize)];
        end
    end
else % nodesize为偶数，取中间2行
    for i=1:node_m
        for j=1:node_n
            nodes{i,j}=img_color((i-1)*nodesize+1:min(i*nodesize,M),(j-1)*nodesize+1:min(j*nodesize,N),:);
            depth(i,j)=mean(img_depth((i-1)*nodesize+1:min(i*nodesize,M),(j-1)*nodesize+1:min(j*nodesize,N)),"all");
            scan_lines{i,j}=[img_depth];
        end
    end
end
% 基本思想：以每个节点对应的空间直线以及节点颜色为标准，进行聚类分析。提取包含底部中间节点的类，生成边界用于导航。
node_last_row=nodes(node_m,:);% 取出最下面一行节点
spatial_line = fittype('a*x + b*y + c', 'coefficients', {'a', 'b', 'c'}, ...
                       'independent', {'x', 'y'}, 'dependent', 'z');
fit_spatial_line = fit(img_depth,spatial_line);