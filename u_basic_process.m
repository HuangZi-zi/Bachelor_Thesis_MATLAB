function output = u_basic_process(img)
% 中值滤波与去除高光

mid_filter_size=9;% 中值滤波器尺寸
min_filter_size=199;% 最小值滤波器尺寸
[h,w,channel]=size(img);

% img=im2double(img);


% 利用最小值滤波得到高光遮罩，并去除高光
R=img(:,:,1);
G=img(:,:,2);
B=img(:,:,3);
min_R = ordfilt2(R,1,ones(min_filter_size,min_filter_size));
min_G = ordfilt2(G,1,ones(min_filter_size,min_filter_size));
min_B = ordfilt2(B,1,ones(min_filter_size,min_filter_size));


ave_R=mean(min_R,"all");
ave_G=mean(min_G,"all");
ave_B=mean(min_B,"all");
% ave=(ave_R+ave_G+ave_B)/3;

min_R(min_R(:)<ave_R)=0;
min_G(min_G(:)<ave_G)=0;
min_B(min_B(:)<ave_B)=0;

mask_R = imerode(R,strel('square',35));
imshow(mask_R);
% 
% % figure();imshow(min_R);
% 
% R1=R-floor(1.*min_R);
% G1=G-floor(1.*min_G);
% B1=B-floor(1.*min_B);
% 
% % R1=R;
% % G1=G;
% % B1=B;
% 
% img1(:,:,1)=R1;
% img1(:,:,2)=G1;
% img1(:,:,3)=B1;
% img1=im2double(img1);
% imshow(im2uint8(img1));
 

% % 用非线性变换解决高光
% img=im2double(img);
% R=img(:,:,1);
% G=img(:,:,2);
% B=img(:,:,3);
% R1 = log(R+1);%灰度图非线性变换，提亮暗部，压暗高光
% G1 = log(G+1);
% B1 = log(B+1);
% 
% img1(:,:,1)=R1;
% img1(:,:,2)=G1;
% img1(:,:,3)=B1;



% 直方图平均
% hist_eq = histeq(gray);
% figure;
% imshow(hist_eq);
% hist_eq_double = im2double(hist_eq);

% 中值滤波
x_dir=img1(:)';
fil_x=medfilt1(x_dir,mid_filter_size);
fil_x_re=pagetranspose(reshape(fil_x,h,w,channel));
y_dir=fil_x_re(:)';
fil_xy=medfilt1(y_dir,mid_filter_size);
fil_xy_re=pagetranspose(reshape(fil_xy,w,h,channel));
img=im2uint8(fil_xy_re);

% x_dir=R1(:)';
% fil_x=medfilt1(x_dir,filter_size);
% fil_x_re=reshape(fil_x,h,w)';
% y_dir=fil_x_re(:)';
% fil_xy=medfilt1(y_dir,filter_size);%2次1维中值0.03024628s
% R1=reshape(fil_xy,w,h)';
% 
% x_dir=G1(:)';
% fil_x=medfilt1(x_dir,filter_size);
% fil_x_re=reshape(fil_x,h,w)';
% y_dir=fil_x_re(:)';
% fil_xy=medfilt1(y_dir,filter_size);%2次1维中值0.03024628s
% G1=reshape(fil_xy,w,h)';
% 
% x_dir=B1(:)';
% fil_x=medfilt1(x_dir,filter_size);
% fil_x_re=reshape(fil_x,h,w)';
% y_dir=fil_x_re(:)';
% fil_xy=medfilt1(y_dir,filter_size);%2次1维中值0.03024628s
% B1=reshape(fil_xy,w,h)';
% 
% img(:,:,1)=R1;
% img(:,:,2)=G1;
% img(:,:,3)=B1;
% 
% img=im2uint8(img);

% 结果可视化
% figure;
% subplot 211;
% imshow(fil_x_re');
% subplot 212;
imshow(img);


% output=fil_xy_re;
output=0;

end