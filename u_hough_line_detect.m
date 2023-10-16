function laneLines = u_hough_line_detect(img, cannyImage)

% Perform Hough Transform to detect lines in the Canny image
[H, theta, rho] = hough(cannyImage);

% Set thresholds for line detection
threshold = 0.5 * max(H(:));
P = houghpeaks(H, 10, 'Threshold', threshold);

% Extract lines from Hough Transform results
lines = houghlines(cannyImage, theta, rho, P, 'FillGap', 20, 'MinLength', 40);

% 绘制所有直线
resultImage = img;
% for i = 1:length(lines)
% xy = [lines(i).point1; lines(i).point2];
% resultImage = insertShape(resultImage, 'Line', xy, 'Color', 'green', 'LineWidth', 1);
% end

for i = 1:length(lines)
    % x=ay+b
    a=(lines(i).point2(1) - lines(i).point1(1)) / (lines(i).point2(2) - lines(i).point1(2));
    b=lines(i).point2(1) - a * lines(i).point2(2);
    lineinfo=struct('a', a, 'b', b, 'point1', lines(i).point1, 'point2', lines(i).point2, 'Merged',false);
    linesinfo(i)=lineinfo;
end

% Define parameters for line merging (adjust as needed)
aThreshold = 0.1; % Maximum slope difference for merging line segments
bThreshold = 20;    % Maximum intercept difference for merging line segments

% Initialize an empty list to store merged lines
mergedLines = [];

% Iterate through each line segment
for i = 1:length(linesinfo)
    currentLine = linesinfo(i);

    % Check if the current line has been merged
    if ~currentLine.Merged
        % Create a new group for merging
        mergeGroup = [currentLine];
        currentLine.Merged = true;

        % Iterate through remaining lines to find merge candidates
        for j = i+1:length(lines)
            candidateLine = linesinfo(j);

            % Check if the candidate line has not been merged and meets merging criteria
            if ~candidateLine.Merged && ...
                    isMergeCandidate(currentLine, candidateLine, aThreshold, bThreshold)
                mergeGroup = [mergeGroup, candidateLine];
                candidateLine.Merged = true;
            end
        end

        % Merge the lines in the group and store the merged line
        mergedLine = mergeLineSegments(mergeGroup);
        mergedLines = [mergedLines; mergedLine];
    end
end

% 绘制合并后所有直线
%     resultImage = img;
%     for i = 1:length(mergedLines)
%     xy = [mergedLines(i).point1; mergedLines(i).point2];
%     resultImage = insertShape(resultImage, 'Line', xy, 'Color', 'red', 'LineWidth', 1);
%     end


% 区分左侧和右侧的线段
leftLaneLines = [];
rightLaneLines = [];
for i = 1:length(mergedLines)
    if mergedLines(i).point1(1) < size(cannyImage, 2) / 2
        if mergedLines(i).a<0
            leftLaneLines = [leftLaneLines; mergedLines(i)];
        end
    else
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
    fprintf("error!!!");
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

    % xx=0:1000;
    % yy=leftLaneLines(sorted_L(1)).a * xx+ leftLaneLines(sorted_L(1)).b;
    %figure;
    imshow(resultImage);
    % plot(xx,yy);
    fprintf("slope of the two lines: %f\t%f\n",leftLaneLines(sorted_L(1)).a,rightLaneLines(sorted_R(1)).a);
    title('直线检测结果');
end


end

%% Function to check if two line segments can be merged
function isMergeable = isMergeCandidate(line1, line2, aThreshold, bThreshold)
% y=ax+b
isMergeable=false;

if abs(line1.a-line2.a)<aThreshold && abs(line1.b-line2.b)<bThreshold
    isMergeable=true;
end
end

%% Function to merge a group of line segments into a single line
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


