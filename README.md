# 使用说明
注：本设计使用Kinect相机，需要安装Kin2驱动（https://github.com/jrterven/Kin2）

## GUI
GUI文件为 GUI.mlapp，打开后正确设置串口com号即可运行。

## 子程序
所有以u_ 开头的文件都是子程序文件。GUI中调用了这些子程序。在不使用GUI的情况下，可以用main_regiongrowing.m进行调用。与GUI实现的功能相同。

u_APF.m：APF路径规划。
u_plane