function output = u_basic_process(img,filtersize,core)
% 中值滤波与去除高光
% core 膨胀和腐蚀的卷积核
mid_filter_size=filtersize;% 中值滤波器尺寸
min_filter_size=filtersize;% 最小值滤波器尺寸
domain=ones(min_filter_size,min_filter_size);
strength=0.21;% 去高光的程度

[h,w,channel]=size(img);

if channel==3 %彩色图，去高光+中值滤波
    % 利用最小值滤波得到高光遮罩，并去除高光
    % 使用HSV空间
    img_hsv=rgb2hsv(img);
    img_v=img_hsv(:,:,3);
    min_fil=ordfilt2(img_v,1,domain);
    mask1=min_fil;
    mask1(mask1(:)<0.8)=0;
    mask1=imgaussfilt(mask1,15);

    mask2=min_fil;
    mask2(mask2(:)<0.95)=0;
    mask2=imdilate(mask2,core);
    mask2=imgaussfilt(mask2,5);

    mask=mask1+mask2;
    img_v=img_v-strength.*mask;
    img_v(img_v<0)=0;
    img_hsv(:,:,3)=img_v;

    img1=hsv2rgb(img_hsv);

    % 均值滤波
    x_dir=img1(:)';
    fil_x=medfilt1(x_dir,mid_filter_size);
    fil_x_re=pagetranspose(reshape(fil_x,h,w,channel));
    y_dir=fil_x_re(:)';
    fil_xy=medfilt1(y_dir,mid_filter_size);
    fil_xy_re=pagetranspose(reshape(fil_xy,w,h,channel));
    img=im2uint8(fil_xy_re);

    output=img;

elseif channel==1 % 深度图，孔洞填充+中值滤波
    img1=im2double(img);
    % 均值滤波
    x_dir=img1(:)';
    fil_x=medfilt1(x_dir,mid_filter_size);
    fil_x_re=pagetranspose(reshape(fil_x,h,w,channel));
    y_dir=fil_x_re(:)';
    fil_xy=medfilt1(y_dir,mid_filter_size);
    fil_xy_re=pagetranspose(reshape(fil_xy,w,h,channel));
    img=im2uint16(fil_xy_re);
    output=img;
end
end