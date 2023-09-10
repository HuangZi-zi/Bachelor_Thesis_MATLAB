clear;
img=imread("demo.png");
imgGray=im2gray(img);
imgBin=imbinarize(img);

processed=u_basic_process(imgGray);
hist=u_histogram(imgGray);
% 
canny=u_find_edge(processed);
% 
 u_find_lane(canny,hist)

