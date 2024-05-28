function laneLines = u_line_hough(img, edge)
resultImage=img;
% 进行hough变换
[H, theta, rho] = hough(edge);

% 设置阈值
threshold = 0.2 * max(H(:));
P = houghpeaks(H, 10, 'Threshold', threshold);

% 提取Hough变换得到的直线
lines = houghlines(edge, theta, rho, P, 'FillGap', 20, 'MinLength', 40);

% 读取直线信息
for i = 1:length(lines)
    % x=ay+b
    a=(lines(i).point2(1) - lines(i).point1(1)) / (lines(i).point2(2) - lines(i).point1(2));
    b=lines(i).point2(1) - a * lines(i).point2(2);
    lineinfo=struct('a', a, 'b', b, 'point1', lines(i).point1, 'point2', lines(i).point2, 'Merged',false);
    linesinfo(i)=lineinfo;
end

% 直线合并的阈值
aThreshold = 0.15;
bThreshold = 50;

mergedLines = [];

% 遍历每个直线段
for i = 1:length(linesinfo)
    currentLine = linesinfo(i);

    % 检查是否合并过
    if ~currentLine.Merged
        mergeGroup = [currentLine];
        currentLine.Merged = true;

        for j = i+1:length(lines)
            candidateLine = linesinfo(j);

            % 检查是否合并过、是否可合并
            if ~candidateLine.Merged && ...
                    isMergeCandidate(currentLine, candidateLine, aThreshold, bThreshold)
                mergeGroup = [mergeGroup, candidateLine];
                candidateLine.Merged = true;
            end
        end

        % 合并直线
        mergedLine = mergeLineSegments(mergeGroup);
        mergedLines = [mergedLines; mergedLine];
    end
end

% 区分左侧和右侧的线段
leftLaneLines = [];
rightLaneLines = [];
for i = 1:length(mergedLines)
    if mergedLines(i).point1(1) < 3*size(edge, 2) / 4 
        if mergedLines(i).a<0
            leftLaneLines = [leftLaneLines; mergedLines(i)];
        end
    elseif mergedLines(i).point2(1) > size(edge, 2) / 4
        if mergedLines(i).a>0 && mergedLines(i).a < Inf
            rightLaneLines = [rightLaneLines; mergedLines(i)];
        end
    end
end

% 根据长度进行排序
lineLengths_L = zeros(length(leftLaneLines), 1);
for i = 1:length(leftLaneLines)
    lineLengths_L(i) = norm(leftLaneLines(i).point1 - leftLaneLines(i).point2);
end
[~,sorted_L] = sort(lineLengths_L, 'descend');

lineLengths_R = zeros(length(rightLaneLines), 1);
for i = 1:length(rightLaneLines)
    lineLengths_R(i) = norm(rightLaneLines(i).point1 - rightLaneLines(i).point2);
end
[~,sorted_R] = sort(lineLengths_R, 'descend');

if isempty(leftLaneLines) && isempty(rightLaneLines)
    fprintf("error!!!\n");
    return;

elseif isempty(leftLaneLines)
    laneLines = [rightLaneLines(sorted_R(1))];
    xy = [rightLaneLines(sorted_R(1)).point1; rightLaneLines(sorted_R(1)).point2];
    resultImage = insertShape(resultImage, 'Line', xy, 'Color', 'cyan', 'LineWidth', 5);
    figure;
    imshow(resultImage);
    title('直线检测结果');
    return;

elseif isempty(rightLaneLines)
    laneLines = [leftLaneLines(sorted_L(1))];
    xy = [leftLaneLines(sorted_L(1)).point1; leftLaneLines(sorted_L(1)).point2];
    resultImage = insertShape(resultImage, 'Line', xy, 'Color', 'cyan', 'LineWidth', 5);
    figure;
    imshow(resultImage);
    title('直线检测结果');
    return;

else
    laneLines = [leftLaneLines(sorted_L(1)),rightLaneLines(sorted_R(1))];
    xy = [leftLaneLines(sorted_L(1)).point1; leftLaneLines(sorted_L(1)).point2];
    resultImage = insertShape(resultImage, 'Line', xy, 'Color', 'cyan', 'LineWidth', 5);
    xy = [rightLaneLines(sorted_R(1)).point1; rightLaneLines(sorted_R(1)).point2];
    resultImage = insertShape(resultImage, 'Line', xy, 'Color', 'cyan', 'LineWidth', 5);

    figure;
    imshow(resultImage);
   
end


end

%% 判断是否可合并的子程序
function isMergeable = isMergeCandidate(line1, line2, aThreshold, bThreshold)
% y=ax+b
isMergeable=false;

if abs(line1.a-line2.a)<aThreshold && abs(line1.b-line2.b)<bThreshold
    isMergeable=true;
end
end

%% 合并直线段的子程序
function mergedLine = mergeLineSegments(lineGroup)
% lineinfo=struct('a', a, 'b', b, 'point1', lines(i).point1, 'point2', lines(i).point2, 'Merged',false);
suma=0;
sumb=0;
l=length(lineGroup);

for i = 1:l
    suma=suma+lineGroup(i).a;
    sumb=sumb+lineGroup(i).b;
    point = struct('p', lineGroup(i).point1,'x', lineGroup(i).point1(1));
    points(2*i-1) = point;
    point = struct('p', lineGroup(i).point2,'x', lineGroup(i).point2(1));
    points(2*i) = point;
end
[~,sorted] = sort([points.x], 'descend');
pointa=points(sorted(1)).p;
pointb=points(sorted(2*l)).p;
mergedLine=struct('a', suma/l, 'b', sumb/l, 'point1', pointa, 'point2', pointb, 'Merged',false);
end
