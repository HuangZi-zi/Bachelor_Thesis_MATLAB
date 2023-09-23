img=imread("Resource\corridor.jpg");
gray=im2gray(img);
BW=imbinarize(gray,"adaptive");

figure
imshow(BW)

