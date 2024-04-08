%% 基于邻域生长的循线
clear;
addpath('C:\Users\YawnFun\Documents\project\Kin2\Mex');
clear
close all

% Create Kinect 2 object and initiqalize it
k2 = Kin2('color','depth');
% i=20;
% 图像尺寸
depth_width = 512; depth_height = 424; outOfRange = 4096;
color_width = 1920; color_height = 1080;

% Create matrices for the images
depth = zeros(depth_height,depth_width,'uint16');
color = zeros(color_height,color_width,3,'uint8');

edges_last_time=zeros(11,3);
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

% 节点大小与确认裁剪的尺寸问题
nodesize=9;
height_cd=375;
width_cd=512;
if rem(height_cd,nodesize)% 如果不能整除，则进一步裁剪
    height_cd=nodesize*floor(height_cd/nodesize);
end
if rem(width_cd,nodesize)
    width_cd = nodesize*floor(width_cd/nodesize)-1;
end
out = zeros(height_cd,width_cd,3,'uint8');
figure, h1 = imshow(out,[]);
title('output');

% 膨胀和腐蚀的卷积核
core=strel('disk',7);

% PID控制器
dir_last_time=0;% 滤波后的，上一次循环的，机器人方向
dir2_last_time=0;% 滤波前的，上一次循环的，预期的下一次，机器人方向
weigh_last_time=0;
kp=1.1;
ki=0.2;
kd=0.1;
integrator=0;

% 语音提示
[audio_stop,Fs]=audioread('Resource\blocked_by_barrier.mp3');


% 在没有运行指令时停止等待
while(size(obj1.UserData,2)~=6) %没有收到指令
    pause(0.1);%等待
end

while(1)
    %tic
    % 检查命令
    if(size(obj1.UserData,2)==6)
        data=obj1.UserData;
        data_s=char(data);
        data_n=str2double(data_s);
        comm=floor(data_n/10000);%命令码
        num=rem(data_n,10000);%命令数
        switch comm
            case 0 % 停止
                sendcomm_stop(obj2);
                fprintf("Stopped by manual input!\n");
                run_flag=0;
            case 1 % 左转
                sendcomm_spin(obj2,comm,num);
            case 2 % 右转
                sendcomm_spin(obj2,comm,num);
            case 3 % 到达
                run_pos=num;
                fprintf("Now arrived at: %d\n",run_pos);
            case 4 % 出发
                run_des=num;
                fprintf("Set destination as: %d\n",run_des);
                run_flag=1;
            case 5 % 调速
                v=num;
                fprintf("New velocity is: %d\n",v);
            otherwise
                fprintf("Wrong Command!\n");
        end
        obj1.UserData=[];
    elseif size(obj1.UserData,2)>6
        fprintf("Wrong Communication!\n");
        obj1.UserData=[];
    end

    %根据命令运行
    if run_flag % 正在运行
        if run_pos==run_des % 到达目的地
            fprintf("Reached destination!\n");
            sendcomm_stop(obj2);
            run_flag=0;
        else % 没有达到
            % 读取新的图像
            validData = k2.updateData;
            if validData

                % Copy data to Matlab matrices
                depth = k2.getDepth;
                color = k2.getColor;
                %         [depth,~]=undistortImage(depth,d_int);
                %         [color,~]=undistortImage(color,c_int);
                color=imresize(color,[375,667]);

                depthColor_c=fliplr(imcrop(color,[89 1 width_cd-1 height_cd-1]));
                depthColor_c=u_basic_process(depthColor_c,9,core);
                depthColor_d=fliplr(imcrop(depth,[1 8 width_cd-1 height_cd-1]));
                depthColor_d(depthColor_d>4096) = 4096;
%                 depthColor_d=depthColor_d.*16;
                depthColor_d=u_basic_process(depthColor_d,9,core);
                %imshowpair(depthColor_d,depthColor_c);
                

                % 查找道路标线
                
                [edges,barrier,barrier_pos]=u_plane_regiongrowing(depthColor_c,depthColor_d,nodesize,core);
                edges=(edges+edges_last_time)/2;% FIR滤波器平滑edges
                if isempty(barrier)% 没有障碍物
                    [out,dir1,dir2,weigh]=u_APF(depthColor_c,edges,v);
                    set(h1,'CData',out);drawnow;
                elseif barrier_pos(1)<width_cd/2 && barrier_pos(2)>width_cd/2 % 停车
                    sendcomm_stop(obj2);
                    sound(audio_stop,Fs);
                    pause(1.5);
                    continue;
                else %避障
                    edges_and_barrier=union(edges,barrier,"rows");
                    [out,dir1,dir2,weigh]=u_APF(depthColor_c,edges_and_barrier,v);
                    set(h1,'CData',out);drawnow;
                end
                edges_last_time=edges;% 更新上一次的边界值
%                 figure(2);imshow(depthColor_d*16);

                %控制器
                new_weigh=weigh/(weigh+weigh_last_time);
                new_weigh_last_time=weigh_last_time/(weigh+weigh_last_time);
                dir_this_time=dir1*new_weigh+dir2_last_time*new_weigh_last_time;
                integrator=integrator+dir_this_time;
                if integrator>v
                    integrator=v;
                elseif integrator<-v
                    integrator=-v;
                end
                dv=fix(kp*dir_this_time+ki*integrator+kd*(dir_this_time-dir_last_time));
                dir_last_time=dir1;
                dir2_last_time=dir2;
                weigh_last_time=weigh;
                
%                 dv=fix(-dir*100);
               sendcomm_run(obj2,v,dv)

            end
        end
    else % 停止运行
        pause(0.1);
    end

%toc
end


% Close kinect object
k2.delete;

close all;
%% 运行函数
function sendcomm_run(port,v,dv)
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
% fprintf("dv:%d\n",dv)
%     fprintf("v:%d \n",v)
vlhex=dec2hex(vl,4);
vrhex=dec2hex(vr,4);
vlh= vlhex(1:2) ;%高位
vll=vlhex(3:4);%低位
vrh= vrhex(1:2) ;%高位
vrl=vrhex(3:4) ;%低位

sendbuff=zeros(1,9);
% 16进制 55 aa 71 04 10 vlh vll vrh vrl sum
sendbuff(1)= 85;
sendbuff(2)= 170;
sendbuff(3)= 113;
sendbuff(4)= 4;
sendbuff(5)= 16;

%             sendbuff(1)= hex2dec('55');
%             sendbuff(2)= hex2dec('aa');
%             sendbuff(3)= hex2dec('71');
%             sendbuff(4)= hex2dec('04');
%             sendbuff(5)= hex2dec('10');
sendbuff(6)= hex2dec(vlh);
sendbuff(7)= hex2dec(vll);
sendbuff(8)= hex2dec(vrh);
sendbuff(9)= hex2dec(vrl);
%校验和
%校验位=前面所有数据之和，取最后两位
add=sum(sendbuff,[1 2 3 4 5 6 7 8 9]);
two_bits_d=rem(add,256);%10进制下对256取余，在16进制下为2位
% two_bits_h= dec2hex(two_bits_d);% 发送数据以10进制存储，因此不需转换
sendbuff(10)= two_bits_d;
write(port,sendbuff,"uint8");
end

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
