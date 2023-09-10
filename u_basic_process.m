function [img] = u_basic_process(inputIMG)
%必须输入灰度图

img=inputIMG;

kernel=strel('square',5);
f_blur=fspecial("gaussian");
%膨胀
img=imdilate(img,kernel);
%模糊
img=imfilter(img,f_blur);
%腐蚀
img=imerode(img,kernel);

% figure;
% imshow(img);
% title("after basic procession")
end