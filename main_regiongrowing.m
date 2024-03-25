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
color = zeros(color_height,color_width,3,'uint8');
out = zeros(375,512,3,'uint8');

k=[];% 键盘输入控制
figure, h = imshow(out,[]);
title('Color2DephAligned (press q to exit)');
set(gcf,'keypress','k=get(gcf,''currentchar'');'); % listen keypress

% calib = k2.getDepthIntrinsics;
% d_int = cameraIntrinsics([calib.FocalLengthX,calib.FocalLengthY], ...
%                          [calib.PrincipalPointX,calib.PrincipalPointY], ...
%                          [depth_height,depth_width]);% 深度相机内参
% k_int = [calib.RadialDistortionSecondOrder,calib.RadialDistortionFourthOrder,calib.RadialDistortionSixthOrder,0,0];
% clear calib
% calib = k2.getColorCalib;
% c_int = cameraIntrinsics([calib.FocalLengthX,calib.FocalLengthY], ...
%                          [calib.PrincipalPointX,calib.PrincipalPointY], ...
%                          [color_height,color_width]);% 彩色相机内参

% 接收读码器信号的端口
clear obj1
obj1=serialport("COM3",9600);
u_QR_Serial(obj1);

% 发送控制指令的端口
clear obj2
obj2=serialport("COM7",19200,'Timeout', 0.2);

% 运行控制
run_flag=0;
run_pos=[];
run_des=[];
v=200;
% 在没有运行指令时停止等待
while(size(obj1.UserData,2)~=6) %没有收到指令
    pause(0.1);%等待
end

while(1)
    if(size(obj1.UserData,2)==6)
        data=obj1.UserData;
        data_s=char(data);
        data_n=str2double(data_s);
        comm=floor(data_n/10000);%命令码
        num=rem(data_n,10000);%命令数
        switch comm
            case 0 % 停止
                sendcomm_stop(obj2);
                break;
            case 1 % 左转
                sendcomm_spin(obj2,comm,num);
            case 2 % 右转
                sendcomm_spin(obj2,comm,num);
            case 3 % 到达
                run_pos=num;
                run_flag=1;
            case 4 % 出发
                run_des=num;
                run_flag=1;
            case 5 % 调速
                v=num;
                run_flag=1;
            otherwise
                fprintf("Wrong Command!\n");
        end
        obj1.UserData=[];
        if (run_flag==1)
            sendcomm_run(k2,obj2,v,run_des,run_pos);
        end
    elseif size(obj1.UserData,2)>6
        fprintf("Wrong Communication!\n");
        obj1.UserData=[];
    else
        sendcomm_run(k2,obj2,v,run_des,run_pos);
    end
end


% Close kinect object
k2.delete;

close all;


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

%% 运行函数
function sendcomm_run(k2,port,v,dst,pos)
if(dst==pos)
    sendcomm_stop;
else
    validData = k2.updateData;
    
    % Before processing the data, we need to make sure that a valid
    % frame was acquired.
    if validData
        % Copy data to Matlab matrices
        depth = k2.getDepth;
        color = k2.getColor; 
%         [depth,~]=undistortImage(depth,d_int);
%         [color,~]=undistortImage(color,c_int);
        color=imresize(color,[375,667]);
        
        depthColor_c=fliplr(imcrop(color,[89 1 511 375]));
        depthColor_c=u_basic_process(depthColor_c);
        depthColor_d=fliplr(imcrop(depth,[1 8 511 374]));
        %imshowpair(depthColor_d,depthColor_c);
        % 
        edges=u_plane_regiongrowing(depthColor_c,depthColor_d);
        [out,dir]=u_APF(depthColor_c,edges);
        figure(1);imshow(out);
%         set(h,'CData',out);
    
%     if ~isempty(k)
%         if strcmp(k,'q')
%             return
%         end
%         
%     end
   


    dv=fix(-dir*100);
    % 增加限幅
    if dv>v
        dv=v;
    elseif dv<-v
            dv=-v;
    end
    vl=v+dv+32768;
    vr=v-dv;
    % 前进为左边反转，右边正转
    % fprintf("left:%d; right:%d\n",vl,vr)
    fprintf("dv:%d\n",dv)
%     fprintf("v:%d \n",v)
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
    write(port,sendbuff,"uint8");


    %     %关闭串口
    %     delete(obj1);
    %     clear obj1;
    %     delete(obj2);
    %     clear obj2;
    end
end
end