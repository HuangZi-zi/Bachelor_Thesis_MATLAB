function numInliers = u_RANSAC(lines)
    % y=mx+b
    % lines: Matrix containing line segments' endpoints (Nx2, where N is the number of segments)
    % inlierThreshold: Maximum distance allowed for an endpoint to be considered an inlier

    % RANSAC（Random Sample Consensus）算法合并线段
    % Define RANSAC parameters
    numIterations = 1000;
    inlierThreshold = 5; % Adjust this threshold as needed
    minNumInliers = 10;  % Minimum number of inliers to consider a line
    
    mergedLines = [];

    for iteration = 1:numIterations
        % Randomly select a line segment
        sample_x = randperm(size(cannyImage,2),2);
        sample_y = randperm(size(cannyImage,1),2);
    
        m = (sample_y(2) - sample_y(1)) / (sample_x(2) - sample_x(1));
        b = sample_y(1) - m * sample_x(1);
    
        % Calculate the coefficients of the line equation (Ax + By + C = 0)
        A = m;
        B = -1;
        C = b;
    
        % Calculate the denominator for distance computation
        denominator = sqrt(A^2 + B^2);
    
        % Initialize the count of inliers
        numInliers = 0;
    
        % Iterate through each line segment
        for i = 1:size(lines, 1)
            % Extract endpoint coordinates
            x = lines(i, 1);
            y = lines(i, 2);
    
            % Calculate the perpendicular distance from the endpoint to the line
            distance = abs(A * x + B * y + C) / denominator;
    
            % Check if the distance is within the inlier threshold
            if distance <= inlierThreshold
                numInliers = numInliers + 1;
            end
        end
    
        % If enough inliers are found, consider it a merged line
        if numInliers >= minNumInliers
            mergedLines = [mergedLines; fittedLine];
        end
    end
    
        % Visualize the merged lines
        imshow(originalImage);
        hold on;
        for i = 1:length(mergedLines)
            plot(mergedLines(i).X, mergedLines(i).Y, 'LineWidth', 2, 'Color', 'r');
        end
        hold off;

end