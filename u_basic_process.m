function output = u_basic_process(img)
% 中值滤波与直方图平均

[h,w,channel]=size(img);

% 转化为灰度图处理
if channel == 3
    gray=rgb2gray(img);
else
    gray = img;
end

gray_double = im2double(gray);

% n1=imnoise(gray_double,'salt & pepper',0.05);
% figure;
% imshow(n1);

gray_log = log(gray_double+1);%灰度图非线性变换，提亮暗部，压暗高光

% 直方图平均
% hist_eq = histeq(gray);
% figure;
% imshow(hist_eq);

% 中值滤波
filter_size=9;
% hist_eq_double = im2double(hist_eq);
% x_dir=hist_eq_double(:)';
x_dir=gray_log(:)';
fil_x=medfilt1(x_dir,filter_size);
fil_x_re=reshape(fil_x,h,w)';
y_dir=fil_x_re(:)';
fil_xy=medfilt1(y_dir,filter_size);%2次1维中值0.03024628s
fil_xy_re=reshape(fil_xy,w,h)';
% filted=medfilt2(n1,[filter_size,filter_size]);%1次2维中值0.18264729s

figure;
subplot 211;
imshow(fil_x_re');
subplot 212;
imshow(fil_xy_re);


output=fil_xy_re;


end