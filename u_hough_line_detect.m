function laneLines = u_hough_line_detect(img, cannyImage)
    % Perform Hough Transform to detect lines in the Canny image
    [H, theta, rho] = hough(cannyImage);
    
    % Set thresholds for line detection
    threshold = 0.5 * max(H(:));
    P = houghpeaks(H, 10, 'Threshold', threshold);
    
    % Extract lines from Hough Transform results
    lines = houghlines(cannyImage, theta, rho, P, 'FillGap', 5, 'MinLength', 100);
    
    % Initialize variables to store lane lines
    leftLaneLines = [];
    rightLaneLines = [];
    
    % Define the angle threshold for classifying lines as left or right lanes
    angleThreshold = 70; % You may need to adjust this value
    
    % Loop through detected lines and classify them as left or right lanes
    for i = 1:length(lines)
        angle = abs(lines(i).theta);
        
        % Filter lines based on their angles
        if angle <= angleThreshold || angle >= (180 - angleThreshold)
            % Line is nearly horizontal, possibly a lane line
            if lines(i).point1(1) < size(cannyImage, 2) / 2
                leftLaneLines = [leftLaneLines; lines(i)];
            else
                rightLaneLines = [rightLaneLines; lines(i)];
            end
        end
    end
    
    % Store the detected lane lines in a structure
    laneLines.left = leftLaneLines;
    laneLines.right = rightLaneLines;

    % Draw the detected lane lines on a copy of the original image
    resultImage = img;
    for i = 1:length( laneLines.left)
        xy = [laneLines.left(i).point1; laneLines.left(i).point2];
        resultImage = insertShape(resultImage, 'Line', xy, 'Color', 'red', 'LineWidth', 2);
    end
    for i = 1:length( laneLines.right)
        xy = [laneLines.right(i).point1; laneLines.right(i).point2];
        resultImage = insertShape(resultImage, 'Line', xy, 'Color', 'red', 'LineWidth', 2);
    end
    
figure;
imshow(resultImage);
title('直线检测结果');
end
