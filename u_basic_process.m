function output = u_basic_process(img,filtersize,core)
% 中值滤波与去除高光

mid_filter_size=filtersize;% 中值滤波器尺寸
min_filter_size=filtersize;% 最小值滤波器尺寸
domain=ones(min_filter_size,min_filter_size);
strength=0.21;% 去高光的程度
% core=strel('disk',7);% 膨胀和腐蚀的卷积核

[h,w,channel]=size(img);

% img=im2double(img);

if channel==3 %彩色图，去高光+中值滤波
    % 利用最小值滤波得到高光遮罩，并去除高光
    

    % --------------------使用HSV空间
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
    % --------------------使用HSV空间---------------------------------END

    % --------------------使用RGB空间
%     gray_img=im2gray(img);
%     R=img(:,:,1);
%     G=img(:,:,2);
%     B=img(:,:,3);
    % min_R = ordfilt2(R,1,ones(min_filter_size,min_filter_size));
    % min_G = ordfilt2(G,1,ones(min_filter_size,min_filter_size));
    % min_B = ordfilt2(B,1,ones(min_filter_size,min_filter_size));
%     min_fil=ordfilt2(gray_img,1,ones(min_filter_size,min_filter_size));

    % ave_R=mean(min_R,"all");
    % ave_G=mean(min_G,"all");
    % ave_B=mean(min_B,"all");
    % ave=(ave_R+ave_G+ave_B)/3;

    % min_R(min_R(:)<200)=0;
    % min_G(min_G(:)<200)=0;
    % min_B(min_B(:)<200)=0;

    % mask_R = imdilate(min_R,strel('disk',5));
    % mask_G = imdilate(min_G,strel('disk',5));
    % mask_B = imdilate(min_B,strel('disk',5));
    % mask_R(mask_R(:)<200)=0;
    % mask_G(mask_G(:)<200)=0;
    % mask_B(mask_B(:)<200)=0;
%     mask=imdilate(min_fil,core);
    % mask=imerode(min_fil,strel('disk',7));
%     mask1=min_fil;
%     mask1(mask1(:)<210)=0;
%     mask1=imdilate(mask1,core);
%     mask1=imgaussfilt(mask1,50);
% 
%     mask2=min_fil;
%     mask2(mask2(:)<230)=0;
%     mask2=imdilate(mask2,core);
%     mask2=imgaussfilt(mask2,10);
% 
%     mask=mask1+mask2;
    % figure();imshow(mask);
    % R1=R-0.2.*mask_R;
    % G1=R-0.2.*mask_G;
    % B1=R-0.2.*mask_B;
    % gray_img1=gray_img-0.2.*mask;

    % % figure();imshow(min_R);
    %
%     R1=R-floor(strength.*mask);
%     G1=G-floor(strength.*mask);
%     B1=B-floor(strength.*mask);
% 
%     img1(:,:,1)=R1;
%     img1(:,:,2)=G1;
%     img1(:,:,3)=B1;
%     img1=im2double(img1);
%     % imshow(im2uint8(img1));

    % --------------------使用RGB空间---------------------------------END

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