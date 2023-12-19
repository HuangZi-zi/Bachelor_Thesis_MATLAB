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
增加了分水岭算法用于分割平面特征。目前效果不理想。验证了ORB特征点检测与匹配的代码，发现难以识别到理想拐角的特征点。修复了图像基本处理函数中的错误，现在真正完成了x-y两次一维中值滤波。拍摄了用于特征提取的照片（shape）

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
更改了层次聚类部分的编程。使用了matlab自带的ahc功能，大幅提升运算速度。修改聚类指标为法向量夹角、空间距离、深度平均值。目前没有做加权。