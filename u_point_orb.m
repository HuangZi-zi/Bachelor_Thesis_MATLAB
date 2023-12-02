img=imread("Resource\conor.png");

blueChannel = img(:, :, 3);
img =im2gray(blueChannel);
img = imgaussfilt(img,2);

% detector=vision.FeatureDetector('ORB');

% keypoints=step(detector, img);

keypoints=detectORBFeatures(img);

imshow(img);
hold on;
plot(keypoints);