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