function [Lane_L_X, Lane_R_X, Lane_Y] = u_find_lane(inputIMG)
% 滑动窗口寻找标线

%转为二值化图像处理

%binaryImage=imbinarize(inputIMG);
binaryImage=inputIMG;

mid=ceil(size(binaryImage, 2)/2);
% [~, left_x_base] = max(hist(1:mid));
% [~, right_x_base] = max(hist(mid+1:end));
% right_x_base=right_x_base+mid;
% midpoint=ceil((left_x_base+right_x_base)/2);

% 定义窗口的高度和数量
windowHeight = size(binaryImage, 1) / 5; % 图像高度除以5
numWindows = 5;
% 定义滑动窗口的宽度和步长
windowWidth = 100; % 窗口宽度（根据曲线宽度选择）
stepSize = 5;     % 滑动步长

% 中心线的坐标
Lane_L_X = [];
Lane_Y = [];
Lane_R_X = [];

% 可视化结果
%figure;
%imshow(inputIMG);
%title('中心线确定结果');
%hold on;

% 滑动窗口分析
for windowIndex = 1:numWindows
    % 计算当前窗口的垂直范围
    startY = round((windowIndex - 1) * windowHeight) + 1;
    endY = round(windowIndex * windowHeight);
    
    % 初始化变量以存储当前窗口内的中心线坐标
    windowvalue=[];
    windowposition=[];

    % 左边
    for x = 1:stepSize:mid
        % 定义窗口的水平范围
        windowX = x:min(x + windowWidth - 1, mid);
        
        % 提取当前窗口的图像
        windowImage = binaryImage(startY:endY, windowX);
        
        % 计算窗口内的平均亮度
        windowBrightness = mean(windowImage(:));
        
        %记录滑动过程中窗口亮度的变化
        windowvalue=[windowvalue,windowBrightness];
        windowposition=[windowposition,round(mean(windowX))];% 中心点 X 坐标
    end
    %寻找窗口与曲线相交最多时
    [~, centerIndex] = max(windowvalue);
    centerX=windowposition(centerIndex);
    centerY = round(mean([startY, endY])); % 中心点 Y 坐标
    % 将当前窗口内的中心线坐标添加到整体中
    Lane_L_X = [Lane_L_X, centerX];
    Lane_Y = [Lane_Y, centerY];

    % 绘制窗口的边界框
    %rectangle('Position', [centerX-windowWidth/2, startY, windowWidth, windowHeight], 'EdgeColor', 'g', 'LineWidth', 5);
    
    % 恢复初始化变量
    windowvalue=[];
    windowposition=[];

    % 右边
    for x = mid:stepSize:size(binaryImage, 2)- windowWidth
        % 定义窗口的水平范围
        windowX = x:min(x + windowWidth - 1, size(binaryImage, 2));
        
        % 提取当前窗口的图像
        windowImage = binaryImage(startY:endY, windowX);
        
        % 计算窗口内的平均亮度
        windowBrightness = mean(windowImage(:));
        
        %记录滑动过程中窗口亮度的变化
        windowvalue=[windowvalue,windowBrightness];
        windowposition=[windowposition,round(mean(windowX))];% 中心点 X 坐标
    end
    %寻找窗口与曲线相交最多时
    [~, centerIndex] = max(windowvalue);
    centerX=windowposition(centerIndex);
%     centerY = round(mean([startY, endY])); % 中心点 Y 坐标
    % 将当前窗口内的中心线坐标添加到整体中
    Lane_R_X = [Lane_R_X, centerX];
%     Lane_R_Y = [Lane_R_Y, centerY];

    % 绘制窗口的边界框
    %rectangle('Position', [centerX-windowWidth/2, startY, windowWidth, windowHeight], 'EdgeColor', 'g', 'LineWidth', 5);
end
    %scatter(Lane_L_X, Lane_Y,'red','filled');
    %scatter(Lane_R_X, Lane_Y,'red','filled');



%hold off;

end