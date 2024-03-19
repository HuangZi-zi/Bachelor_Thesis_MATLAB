%% 基于邻域生长的循线
clear;
addpath('C:\Users\YawnFun\Documents\project\Kin2\Mex');
clear
close all

% Create Kinect 2 object and initiqalize it
% Available sources: 'color', 'depth', 'infrared', 'body_index', 'body',
% 'face' and 'HDface'
k2 = Kin2('color','depth');
i=20;
% images sizes
depth_width = 512; depth_height = 424; outOfRange = 4096;
color_width = 1920; color_height = 1080;

% Create matrices for the images
depth = zeros(depth_height,depth_width,'uint16');
depthColor = zeros(depth_height,depth_width,3,'uint8');

k=[];% 键盘输入控制
figure, h = imshow(depthColor,[]);
title('Color2DephAligned (press q to exit)');
set(gcf,'keypress','k=get(gcf,''currentchar'');'); % listen keypress

calib = k2.getDepthIntrinsics;
d_int = cameraIntrinsics([calib.FocalLengthX,calib.FocalLengthY], ...
                         [calib.PrincipalPointX,calib.PrincipalPointY], ...
                         [depth_height,depth_width]);% 深度相机内参
k_int = [calib.RadialDistortionSecondOrder,calib.RadialDistortionFourthOrder,calib.RadialDistortionSixthOrder,0,0];
calib = k2.getColorCalib;
c_int = cameraIntrinsics([calib.FocalLengthX,calib.FocalLengthY], ...
                         [calib.PrincipalPointX,calib.PrincipalPointY], ...
                         [color_height,color_width]);% 彩色相机内参

while true
    % Get frames from Kinect and save them on underlying buffer
    validData = k2.updateData;
    
    % Before processing the data, we need to make sure that a valid
    % frame was acquired.
    if validData
        % Copy data to Matlab matrices
        depth = k2.getDepth;
        color = k2.getColor; 
        [depth,~]=undistortImage(depth,d_int);
        [color,~]=undistortImage(color,c_int);
        color=imresize(color,[375,667]);
        depthColor_d=imcrop(depth,[1 8 374 511]);
        depthColor_c=imcrop(color,[1 89 374 511]);
        imshowpair(depthColor_d,depthColor_c);
        %set(h,'CData',depthColor);
        %u_plane_regiongrowing(depthColor,depth);

    end
    if ~isempty(k)
        if strcmp(k,'q')
            break
        end
        
    end
    pause(0.02)
end

% Close kinect object
k2.delete;

close all;

        % img=imread("Resource\mainline.jpg");
        % img=imread("Resource\turn.jpg");
        v=200;
        fismat=readfis('yz001');
        % cam=webcam(2);

        % 接收读码器信号的端口
%         clear obj1
%         obj1=serialport("COM10",9600);
%         u_QR_Serial(obj1);

        % 发送控制指令的端口
        %     clear obj2
        %     obj2=serialport("COM9",19200,'Timeout', 0.2);

        for iii=1:500
            % img = snapshot(cam);
            % img=u_segment(img);
            img =imcrop(img,[90,150,1100,400]);
            processed=u_basic_process(img);
            % hist=u_histogram(img);
            edge=u_find_edge(img);
            % [Lane_L_X, Lane_R_X, Lane_Y]=u_slide_window(canny,hist);
            % u_fit(Lane_L_X, Lane_R_X, Lane_Y, img);

            laneLines=u_line_hough(img, edge);

            % 计算最下面三行的中心座标
            c3=0;
            for i=399:401
                % x=ay+b
                c3=c3+(i*laneLines(2).a+laneLines(2).b) + (i*laneLines(1).a+laneLines(1).b);
            end
            c3=c3/6;
            fprintf("c3:%d\n",c3)
            output= evalfis(fismat,c3);% 模糊控制
            % 确定输入为位置偏差，其论域范围确定为[−6,6]，而输出为转速，其论域范围为[−10,10]
            dv=-fix(output*20);
            vl=v+dv+32768;
            vr=v-dv;
            % 前进为左边反转，右边正转5
            % fprintf("left:%d; right:%d\n",vl,vr)
            % fprintf("dv:%d\n",dv)
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
%             write(obj2,sendbuff,"uint8");
        end

        % 停止机器人
        sendcomm_stop(obj2)
        %关闭串口
%         delete(obj2);
%         clear obj2;



%% 停止函数
function sendcomm_stop(port)
    sendbuff=zeros(1,9);
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
    write(port,sendbuff,"uint8")
end

%% 旋转函数
function sendcomm_spin(port,dir,time)
    %    01，左转（地面
    %    02，右转（地面
    v=200;
    vl=v+32768*(dir-1);
    vr=v+32768*(dir-1);
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
    sendbuff(10)= two_bits_d;
    write(port,sendbuff,"uint8");

    pause(time/10);

    % 停止机器人
    sendcomm_stop(port);
end
