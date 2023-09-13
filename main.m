clear;
img=imread("Resource\demo.png");
imgGray=im2gray(img);
imgBin=imbinarize(img);

processed=u_basic_process(imgGray);
hist=u_histogram(imgGray);
canny=u_find_edge(processed);
[Lane_L_X, Lane_R_X, Lane_Y]=u_find_lane(canny,hist);
u_fit(Lane_L_X, Lane_R_X, Lane_Y, img);

imgCode=imread("Resource\QRcode.jpg");
