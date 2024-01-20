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