function output = u_plane_watershed(img)

[h,w,channel]=size(img);
% 转化为灰度图处理
if channel == 3
    gray=rgb2gray(img);
else
    gray = img;
end

gray=u_basic_process(gray);

%求梯度
gmag=imgradient(gray);

%标记前景
% using "opening-by-reconstruction" and "closing-by-reconstruction" to "clean" up the image
se = strel('disk',20);
Ie = imerode(gray,se);
Iobr = imreconstruct(Ie,gray);

Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);%输出需要取反
fgm = imregionalmax(Iobrcbr);%使用局部极大值标记前景
se2 = strel(ones(5,5));
fgm2 = imclose(fgm,se2);
fgm3 = imerode(fgm2,se2);%处理前景标记
fgm4 = bwareaopen(fgm3,80);%去除少于某个阈值的像素块

%标记背景
bw = imbinarize(Iobrcbr);%将前景二值化，黑色的为背景
D = bwdist(bw);%距离变换。像素值=该像素到最近的非零像素的距离
DL = watershed(D);%对距离变换求分水岭
bgm = DL == 0;%取"山脊"，用于将背景与前景分割开

%计算输出
gmag2 = imimposemin(gmag, bgm | fgm4);
L = watershed(gmag2);
Lrgb = label2rgb(L,'jet','w','shuffle');
figure
imshow(gray);
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;
title('Colored Labels Superimposed Transparently on Original Image')
output=0;
end