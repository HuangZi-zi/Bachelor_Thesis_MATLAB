function u_QR_Serial(s)

configureTerminator(s,"LF");% 结束符
flush(s);% 清空数据
configureCallback(s,"byte",1,@readSerialData);
end

function readSerialData(src,~)
data = fread(src,1);
src.UserData(end+1) = data;
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