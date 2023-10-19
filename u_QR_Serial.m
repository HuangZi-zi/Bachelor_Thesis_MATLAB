clear;
s=serialport("COM10",9600);
% pause(5);
% read(s,9,"char")
configureTerminator(s,"LF");
flush(s);
% s.UserData = struct("Data",[],"Count",1);
configureCallback(s,"byte",1024,@readSerialData);
count=0;
for i=1:500
    pause(0.01);
    if(strcmp(s.UserData,'400'))
        fprintf("Yes! %d\n",count)
        s.UserData='';
    end
    count=count+1;
end
% clear s;

function readSerialData(src,~)
    data = fscanf(src);
    src.UserData = data;
%     fprintf(data);
end