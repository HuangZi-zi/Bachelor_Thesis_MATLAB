%%
% clear;
% R=1;
% 
% v_desired=1.1*sqrt((5-0)^2+(5-0)^2)/10;% 经验值，自拟速度
% 
% t1=0.125*pi*2*R/v_desired;% 第一段圆弧
% v_l(1:floor(t1/0.01))=1.2*v_desired;
% v_r(1:floor(t1/0.01))=0.8*v_desired;
% 
% t2=sqrt((5-1-0-0)^2+(5-1-0-0)^2)/v_desired;% 直线
% v_l(floor(t1/0.01)+1:floor(t1/0.01)+1+floor(t2/0.01))=v_desired;
% v_r(floor(t1/0.01)+1:floor(t1/0.01)+1+floor(t2/0.01))=v_desired;
% 
% t3=t1; % 第二段圆弧
% v_l(floor(t1/0.01)+1+floor(t2/0.01)+1:floor(t1/0.01)+1+floor(t2/0.01)+1+floor(t3/0.01))=1.2*v_desired;
% v_r(floor(t1/0.01)+1+floor(t2/0.01)+1:floor(t1/0.01)+1+floor(t2/0.01)+1+floor(t3/0.01))=0.8*v_desired;
% 
% v_l(length(v_l):1000)=0;% 停止
% v_r(length(v_r):1000)=0;
% t=0:0.01:length(v_r)/100-0.01;
% 
% v_l=[t', v_l'];% 第一列为时间戳，第二列为数据
% v_r=[t', v_r'];

%%
% clear;
% W=0.4;
% Ts=0.01;
% % startup_rvc % 启动机器人工具箱
% 
% [x,dx,ddx]=tpoly(0,5,10/Ts,0,0.001);
% [y,dy,ddy]=tpoly(0,5,10/Ts,0.001,0);
% 
% % plot(x,y);
% tanphi=dx./dy;
% phi=atan(tanphi);
% dphi(2:1000)=diff(phi);
% 
% dphi=dphi';
% 
% v_desire=sqrt(dx.^2+dy.^2)/Ts;
% 
% t=0:Ts:10-Ts;
% % v_desire=ones(1,1000);
% 
% v_l=[t',  (v_desire+W.*dphi./Ts./2)];
% v_r=[t',  (v_desire-W.*dphi./Ts./2)];
% 
% % v_r(2,2)=0;
% % v_l(2,2)=0;
% % v_r=(v+W/2*dphi)';
% % v_l=(v-W/2*dphi)';

%% 从simulink读取数据并绘图（跟踪直线）
x_out=out.x_out.data;
y_out=out.y_out.data;
plot(x_out,y_out,"-.");xlabel("x/m");ylabel("y/m");
hold on;
grid on;
a=out.lane.a.data;
b=out.lane.b.data;
c=out.lane.c.data;
point_x=min(x_out):1/size(x_out,1):max(x_out);
point_y=-a.*point_x./b-c;
plot(point_x,point_y,"-");
legend('trajectory','lane');

%% 从simulink读取数据并绘图（跟踪任意曲线）
x_out=out.x_out.data;
y_out=out.y_out.data;
plot(x_out,y_out);
hold on;
grid on;
x_in=out.lane.data(:,1);
y_in=out.lane.data(:,2);
plot(x_in,y_in,".");

legend('trajectory','lane');