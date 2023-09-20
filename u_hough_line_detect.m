function laneLines = u_hough_line_detect(img, cannyImage)
    
    % Perform Hough Transform to detect lines in the Canny image
    [H, theta, rho] = hough(cannyImage);
    
    % Set thresholds for line detection
    threshold = 0.5 * max(H(:));
    P = houghpeaks(H, 10, 'Threshold', threshold);
    
    % Extract lines from Hough Transform results
    lines = houghlines(cannyImage, theta, rho, P, 'FillGap', 5, 'MinLength', 100);
    
    % 区分左侧和右侧的线段
    leftLaneLines = [];
    rightLaneLines = [];
    for i = 1:length(lines)
        if lines(i).point1(1) < size(cannyImage, 2) / 2
            leftLaneLines = [leftLaneLines; lines(i)];
        else
            rightLaneLines = [rightLaneLines; lines(i)];
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

    
    
    % Draw the detected lane lines on a copy of the original image
    resultImage = img;
    xy = [leftLaneLines(sorted_L(1)).point1; leftLaneLines(sorted_L(1)).point2];
    resultImage = insertShape(resultImage, 'Line', xy, 'Color', 'red', 'LineWidth', 2);
    xy = [rightLaneLines(sorted_R(1)).point1; rightLaneLines(sorted_R(1)).point2];
    resultImage = insertShape(resultImage, 'Line', xy, 'Color', 'red', 'LineWidth', 2);
    
    
    % Store the detected lane lines in a structure
    laneLines = [leftLaneLines(sorted_L(1)),rightLaneLines(sorted_R(1))];
    
    figure;
    imshow(resultImage);
    title('直线检测结果');
end
