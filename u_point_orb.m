function output = u_point_orb(img,ideal)

% img=imread("Resource\conor.png");
% imread("Resource\luckin_1.jpg");
img_blueChannel = img(:, :, 3);
img =im2gray(img_blueChannel);
img = imgaussfilt(img,1);
ideal_blueChannel = ideal(:, :, 3);
ideal =im2gray(ideal_blueChannel);
ideal = imgaussfilt(ideal,1);

keypoints1=detectORBFeatures(img);
% newPoints = selectUniform(keypoints,20,size(img));
[features1,valid_points1] = extractFeatures(img,keypoints1);

keypoints2=detectORBFeatures(ideal);
% newPoints2 = selectUniform(keypoints2,20,size(ideal));
[features2,valid_points2] = extractFeatures(ideal,keypoints2);

% figure;
% imshow(img);
% hold on;
% plot(keypoints1);
% figure;
% imshow(ideal);
% hold on;
% plot(keypoints2);

indexPairs = matchFeatures(features2,features1);
matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);

figure; ax = axes;
showMatchedFeatures(ideal,img,matchedPoints2,matchedPoints1,"montag",Parent=ax);
legend(ax,"Matched points 1","Matched points 2");

% 
% output=newPoints;

end