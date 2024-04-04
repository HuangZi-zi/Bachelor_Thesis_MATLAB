# IMGprocess_MATLAB
利用matlab完成图像处理

## 0909
完成对俯视图的中心线提取

## 0910
完成图像处理和边缘识别，目前右侧中心线提取失效

## 0911
解决了右侧中心线提取失效的问题

## 0912
完成了二次函数拟合

## 0913
库文件来自https://linuxtv.org/downloads/zbar/binaries/

## 0919
修改了函数结构，现只需传入原始图像即可

新增使用Hough直线检测的模块，在走廊场景实现了识别

## 0920
新增了道路场景。优化了直线筛选逻辑，在道路场景实现了识别。

## 0921
新增了室内场景（square）

新增了角度过滤规则

## 0922
新增了基于RGB三通道的边缘检测，但在室内场景效果仍不理想

## 0923
新增直线合并，目前缺少一个合并后直线端点的确定算法

## 0924
完成了直线合并算法。使用自适应二值化易造成直线的小幅度偏转失准。考虑再加一层canny边缘检测

## 0926
新增了左右失去标线的判别

## 1010
完成了基于滑动窗口的循线。该方法的问题在于拟合直线参数的微小变化可能导致计算结果的大幅度变化，从而导致不稳定。

## 1016
整理了先前的工程，完成了基于Hough检测的循线编程部分

## 1017
尝试整合二维码识别，发现只能识别有3个定位点的二维码，且容易出现错误。考虑更换其他的识别程序。

## 1019
选择采用点化+聚类分析的方法。目前使用圆形滑窗，运算速度慢。

## 1020
完成串口中断的调试，现在可以利用串口读入信息来调节速度

## 1021
完成测试，在开环条件下了解了机器人转向角度与速度、时间的关系

## 1022
完成二维码控制小车运动程序编写

## 1028
修改控制逻辑，不再直接进入运行的程序，而由二维码输入控制启动

## 1202
开始了ORB特征点提取

## 1205
增加了.gitignore，完成了动态ROI的创立部分，还差动态调整部分

## 1206
增加了分水岭算法用于分割平面特征。目前效果不理想。
验证了ORB特征点检测与匹配的代码，发现难以识别到理想拐角的特征点。
修复了图像基本处理函数中的错误，现在真正完成了x-y两次一维中值滤波。拍摄了用于特征提取的照片（shape）

## 1208
增加了快速平面检测算法。目前完成了部分编程。待解决的问题包括每一个节点拟合直线时的数据传输问题。

## 1209
完成了快速平面检测算法中初始化的部分。待解决的问题包括聚类分析法部分的编程。

## 1210
编写了主成分分析法拟合平面。修改了算法逻辑为非函数形式。优化了图像初始化部分的编程。开始了层次聚类分析的编程。

没有完成：rejectedge部分判断是否处于“墙角”的部分

## 1211
完成了层次聚类部分的编程。目前没有做edge部分的判别。目前程序执行速度很慢。

## 1219
更改了层次聚类部分的编程。使用了matlab自带的ahc功能，大幅提升运算速度。
修改聚类指标为法向量夹角、空间距离、深度平均值。目前没有做加权。

## 1220
更改了存储法向量的数据结构为矩阵，期望可以提高运算速度。使用实拍图片进行测试，目前效果不理想。
主要问题在于图像的精度较低，边缘不明显。

## 0106
新思路：利用深度图的梯度信息进行聚类分析。进行实验，效果不佳。
新思路：从图像底部中心开始生长，遇到一定条件则停止，可以提取出地面，供导航使用。
新增了聚类法平面提取的函数文件，预期生成c++代码供点云数据使用。

## 0108
完成了从底部中心开始生长的边缘检测部分，并通过了测试，效果良好。新增了由边缘检测结果进一步检验平面的程序，还未进行测试。
目前并不是严谨的邻域生长法，因为为了节约运算时间，从底部运算结束后，就只取了两个端点进行下一步运算，
因此可能会出现内部边缘而未检测到的问题。

## 0109
安装了kinect驱动，现在将所有操作全部移植到matlab。
完成了基于深度与色彩对齐图上直线特征的平面校验。下一步问题在于删除无效区域，计算道路中心。
下一步需要比较各种方法，完成创新点的立论。

## 0110
修复了软件bug。测试一幅图像处理时间约为0.0427s。测试了使用点云+聚类分析分割平面，用时约为每帧1s

## 0111
测试了点提取的代码，采集了新的角点检测测试图像，并绘制了模板。发现目前方法的视角不变性较差，考虑使用其他方法做点提取。
定义了基于像素色彩和明度的相似性指标，并对应修改了线提取的代码。

## 0120
新增了用于测试的图片

## 0224
假期工作备份。包括针对彩色图像的孔洞填充和滤波、手撮编程的SIFT等。

## 0225
0225进度汇报新增的 去高光操作的限幅。效果一般。

## 0305
关于新增CAPE算法的matlab移植，因为数组下标出现问题，没有继续。
关于反投影回点云，并利用点云的xyz进行聚类分析。效果不好。

## 0306
修改了hough直线检测中的bug。之前的程序中在筛选左侧标线和右侧标线存在判断条件不准确的问题，现在可以正确地提取出最长的的直线。

## 0307
新增了一维中值滤波与二维中值滤波的对比（在draft中）。修改了邻域生长法绘图的输出部分（在u_plane_regiongrowing中）。

## 0308
新增了曲线场景的边缘线提取，测试效果邻域生长法可用。

## 0314
将机器人仿真的相关文件移动到了本项目中。
新增了人工势场法路径规划算法。编程还在进行中。

## 0318
完成了人工势场法路径规划的相关文件，测试效果较好

## 0319
使用图像空间的裁剪与缩放完成了深度图与彩色图的对齐。由于深度图有孔洞，造成坐标映射时超过彩色图范围，造成读取的对齐结果有大量缺失
使用本方法不存在这个问题，可以获得原始质量的彩色图，且误差仍然可控。

## 0320
完成了基于邻域生长和人工势场法的路径规划主程序，由于硬件原因未测试成功。机器人充电完成后做测试。

## 0321
完成了基本的测试，发现AHF算法有问题，导致路径规划会跑出边界。重构了部分计算，修改了参数，效果得到提升。

## 0324
修改了AHF算法的逻辑，消除了超出边界的问题。开始目前的主要问题在于高光干扰。已经验证非线性映射、最小值滤波遮罩不好用。考虑怎加形态学操作。

## 0325
使用腐蚀得到了较好的高光遮罩，基本可以处理高光问题。目前的问题是边缘检测有时会跑飞。开始了新边缘检测算法的编程。

## 0326
完成了按节点拟合直线的编程。目前在识别障碍物下仍有问题，考虑使用新的原则进行障碍物的处理。

## 0327
完成了障碍物避让的程序，测试后基本符合预期。构思了一种衡量同一行节点间拟合直线相似程度的方法，待验证。
收集了存在障碍物情况下的图片用于测试。

## 0403
修改了二维码通信部分，完成了全部二维码控制的功能（启动，停止，目的地判别，调速等），且程序逻辑更简单。将循线程序放回主函数中，避免大量参数的传递。

## 0404
进一步修改了二维码通信部分。将卷积核的创建部分单独提出来，修改了结果图像的显示方式，将算法运行一次的时间从0.15s左右减少到0.09s左右。
目前的问题为，边缘检测结果不稳定，造成机器人震荡，考虑增加滤波或者pid。