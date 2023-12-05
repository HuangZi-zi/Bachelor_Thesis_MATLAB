function u_ROI(laneLines,img)

% x=ay+b

a1=laneLines(1).a;
a2=laneLines(2).a;

b1=laneLines(1).b;
b2=laneLines(2).b;

width_dyn=40;%动态的ROI宽度

% img=imread("Resource\mainline.jpg");

[h,w,channel]=size(img);

line_start_x=0;
line_start_y1=round(-b1/a1);
line_start_y2=round(-b2/a2);

line_end_x=w;
line_end_y1=round((line_end_x-b1)/a1);
line_end_y2=round((line_end_x-b2)/a2);

xi=[line_start_x,line_start_x,line_end_x,line_end_x];
y1i=[line_start_y1-width_dyn,line_start_y1+width_dyn,line_end_y1+width_dyn,line_end_y1-width_dyn];
y2i=[line_start_y2-width_dyn,line_start_y2+width_dyn,line_end_y2+width_dyn,line_end_y2-width_dyn];

mask1=roipoly(img,xi,y1i);
mask2=roipoly(img,xi,y2i);
mask=bitor(mask1,mask2);
roi = bsxfun(@times, img, cast(mask, class(img)));

figure;
subplot(1, 2, 1);
imshow(img);
title('Original Image');

subplot(1, 2, 2);
imshow(roi);
title('Selected Region');
end