%% 基于滑动窗口的双线循线
if(0)
    clear;
    cam=webcam(2);
    fismat=readfis('yz001');
    v=200;
    % 接收读码器信号的端口
    clear obj1
    obj1=serialport("COM10",9600);
    u_QR_Serial(obj1);

    % 发送控制指令的端口
    clear obj2
    obj2=serialport("COM13",19200,'Timeout', 0.2);
    
    % 运行控制
    run_flag=0;
    run_pos=[];
    run_des=[];
    while(size(obj1.UserData,2)~=6) %没有收到指令
       pause(0.1);%等待
    end

    while(1)
        if(size(obj1.UserData,2)==6)
            data=obj1.UserData;
            data_s=char(data);
            data_n=str2num(data_s);
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
                sendcomm_run(cam,obj2,fismat,v,run_des,run_pos);
            end
        elseif size(obj1.UserData,2)>6
            fprintf("Wrong Communication!\n");
            obj1.UserData=[];
        else
            sendcomm_run(cam,obj2,fismat,v,run_des,run_pos);
        end
    end
end

%% 基于Hough检测的双线循线

if(1)
        clear;
        img=imread("Resource\mainline.jpg");
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
            % 前进为左边反转，右边正转
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
end

%% 滑动窗口循线运行
function sendcomm_run(cam,port,fismat,v,dst,pos)
if(dst==pos)
    sendcomm_stop;
else
    img = snapshot(cam);
    % imshow(i);
    % preview(cam)
    % pict =imcrop(i,[0,size(i,1)*0.25,size(i,2),size(i,1)*0.25]);
    pict =imcrop(img,[90,150,1100,400]);
    im=rgb2gray(pict);% 转为灰度
    im = histeq(im);% 直方图平均
    w=3;
    bk=double(im);
    bk=ordfilt2(bk,1,ones(w,w) ,'symmetric');% 最小值滤波
    h=ones(w,w)/(w*w);
    bk=imfilter(bk,h,'replicate');% 均值滤波
    newIm=imsubtract(double(im),bk);% im- bk（小于0的截断为0）
    th=graythresh(newIm/255)*0.40;
    newIm=im2bw(newIm/255,th);% 二值化
    img_fill = imfill(newIm, 'holes');% 孔洞填充
    img = not(img_fill);% 反色
    c = imopen( img, ones(4,4) );% 开运算
    c= bwareafilt(c, 2);% 提取对象，只保留面积最大的一个
    c=imclose(c, ones(32,32) );% 闭运算
    % d = edge( c ,'sobel');% 边缘检测
    d=c;
    % hist=u_histogram(img);
    % canny=u_find_edge(processed);
    [Lane_L_X, Lane_R_X, Lane_Y]=u_slide_window(d);
    c3=u_fit(Lane_L_X, Lane_R_X, Lane_Y, d);
    fprintf("c3:%d\n",c3)

    output= evalfis(fismat,c3);% 模糊控制
    % 确定输入为位置偏差，其论域范围确定为[−6,6]，而输出为转速，其论域范围为[−10,10]
    dv=-fix(output*20);
    vl=v+dv+32768;
    vr=v-dv;
    % 前进为左边反转，右边正转
    % fprintf("left:%d; right:%d\n",vl,vr)
    % fprintf("dv:%d\n",dv)
    fprintf("v:%d \n",v)
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
