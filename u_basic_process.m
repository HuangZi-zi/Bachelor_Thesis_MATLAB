function outputIMG = u_basic_process(inputIMG)


% 转化为灰度图处理
IMGgray=rgb2gray(inputIMG);

kernel=strel('square',5);

%膨胀
IMGdil=imdilate(IMGgray,kernel);
%模糊
IMGblur = imgaussfilt(IMGdil,8);
%腐蚀
IMGero=imerode(IMGblur,kernel);

outputIMG=IMGero;

figure;
imshow(outputIMG);
title("图像基本处理")

end