clear;
img=imread("Resource\corridor.jpg");
img=u_segment(img);

% processed=u_basic_process(img);
% hist=u_histogram(img);
  canny=u_find_edge(img);
%  [Lane_L_X, Lane_R_X, Lane_Y]=u_find_lane(canny,hist);
%  u_fit(Lane_L_X, Lane_R_X, Lane_Y, img);

laneLines=u_hough_line_detect(img, canny);

