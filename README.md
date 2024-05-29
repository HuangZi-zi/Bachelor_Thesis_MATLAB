# 使用说明
注：本设计使用Kinect相机，需要安装Kin2驱动（https://github.com/jrterven/Kin2）

## 主程序
GUI文件为 GUI.mlapp，打开后正确设置串口com号即可运行。

在不使用GUI的情况下，运行main_regiongrowing.m，与GUI实现的功能相同。

## 子程序
所有以u_ 开头的文件都是子程序文件。主程序中调用了这些子程序。

u_APF.m：APF路径规划。
u_basic_process.m：图像基本处理（均值与去高光）。
u_edge.m：sobel算子提取边缘。
u_line_hough.m：hough变换及优化筛选。
u_plane_regiongrowing.m：邻域生长法线特征提取、利用面特征优化线特征、障碍物提取。
u_QR_Serial.m：二维码识读器串口中断接收配置。

## 仿真程序
所有仿真程序均在simulation文件夹内。
my_robot.slx：AGV跟踪直线、曲线轨迹的仿真环境。
final.m：AGV跟踪直线轨迹的绘图程序。
tune_pid.slx：系统辨识与PID校正仿真环境。

## 其他文件
draw_pictures.m文件为论文写作过程中绘图所用，实际数据和绘图有所不同，仅供参考。
Resource文件夹内保存有测试所用图片和部分测试结果。