%% 从simulink读取数据并绘图（跟踪直线）
x_out=out.x_out1.data;
y_out=out.y_out1.data;
plot(x_out,y_out,"-.");xlabel("x/m");ylabel("y/m");
hold on;
grid on;
a=out.lane1.a.data;
b=out.lane1.b.data;
c=out.lane1.c.data;
point_x=min(x_out):1/size(x_out,1):max(x_out);
point_y=-a.*point_x./b-c;
plot(point_x,point_y,"-");
legend('trajectory','lane');

