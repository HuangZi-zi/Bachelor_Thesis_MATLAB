function colHistograms = u_histogram(inputIMG)
% 转为灰度图处理
inputIMG=im2gray(inputIMG);

% 初始化一个数组来存储每一行的直方图
colHistograms = zeros(size(inputIMG, 1), 1);
% 阈值
whiteThreshold = 200;

% 遍历每一列并计算白色像素数量
for col = 1:size(inputIMG, 2)
    % 使用阈值将图像二值化，只保留白色像素
    binaryCol = inputIMG(:, col) < whiteThreshold;
    colHistograms(col) = sum(binaryCol);
end
% 可视化按列白色像素数量
% plot(colHistograms);
% xlabel('列');
% ylabel('白色像素数量');
% title('按列白色像素数量直方图');
end