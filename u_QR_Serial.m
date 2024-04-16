function u_QR_Serial(s)
% clear;
%%释放串口资源
% % List all available serial ports
% portList = serialportlist;
% % Close each open serial port
% for i = 1:numel(portList)   
%         s = serial(portList{i});
%         fclose(s);
%         delete(s);
% end
% delete(instrfindall)%和这一句的作用相同


configureTerminator(s,"LF");% 结束符
flush(s);% 清空数据
%s.UserData = struct("Data",[],"Count",0);
configureCallback(s,"byte",1,@readSerialData);
% count=0;
% for i=1:500
%     pause(0.01);
%     if(strcmp(s.UserData,'400'))
%         fprintf("Yes! %d\n",count)
%         s.UserData='';
%     end
%     count=count+1;
% end
% clear s;

% 释放串口资源

% % 创建串口资源
% Com_Obj = serial('COM10','BaudRate',9600);
% % 设置串口中断,接收到1个字节中断
% set(Com_Obj,'BytesAvailableFcnMode','byte');
% set(Com_Obj,'BytesAvailableFcnCount', 1);
% % 设置窗口接收回调
% Com_Obj.BytesAvailableFcn=@ComReceiveCallback;
% % 查询串口是否已经打开
% if(Com_Obj.Status == "closed")
%     % 打开串口
%     fopen(Com_Obj);
%     if(Com_Obj.Status == "open")
%         fprintf(1,'串口打开成功\n'); 
%     else
%         fprintf(1,'串口打开失败\n');     
%     end  
% else
%     fprintf(1,'串口被占用\n'); 
% end 
% 
% % 串口接收数据
% function ComReceiveCallback(Com_Obj, src, event)
% Com_DataReceive = fread(Com_Obj,1);
% disp(Com_DataReceive);
% % 用户协议包解析
%             
% end
end

function readSerialData(src,~)
    data = fread(src,1);
    src.UserData(end+1) = data;
    %src.UserData.Count = src.UserData.Count + 1;
%     disp(data);   
end



%% QR Code通信协议：
% 1. 命令类型
%     00，停止（手动输入）
%     01，左转（地面
%     02，右转（地面
%     03，到达（地面
%     04，出发（手动
%     05，调速（手动
% 2. 命令参数
%     当命令为左转或右转时：xxxx，以200为速度旋转的时间，eg：6.5s->65->0065
%     当命令为到达时：xxxx，当前位置代号，eg：A地->0000
%     当命令为出发时：xxxx，下一个要到达的位置代号，eg：B地->0001
%     当命令为调速时：xxxx，基本速度，eg：0400