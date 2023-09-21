function outputIMG = u_basic_process(inputIMG)


% 转化为灰度图处理
IMGgray=rgb2gray(inputIMG);

kernel=strel('square',10);

%模糊
IMGblur = imgaussfilt(IMGgray,4);

%膨胀
IMGdil=imdilate(IMGblur,kernel);

%腐蚀
IMGero=imerode(IMGdil,kernel);

outputIMG=IMGero;

figure;
imshow(outputIMG);
title("图像基本处理")

end