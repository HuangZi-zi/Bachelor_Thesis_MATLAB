%% 输出裁剪后的图片
nodesize=9;
height_cd=375;
width_cd=512;
if rem(height_cd,nodesize)% 如果不能整除，则进一步裁剪
    height_cd=nodesize*floor(height_cd/nodesize);
end
if rem(width_cd,nodesize)
    width_cd = nodesize*floor(width_cd/nodesize)-1;
end
color=imread("Resource\snapc2.png");% curl.jpg %snapc30.png %
color=imresize(color,[375,667]);
% depthColor_c=fliplr(imcrop(color,[89 1 width_cd-1 height_cd-1]));
depthColor_c=imcrop(color,[89 1 width_cd-1 height_cd-1]);
imshow(depthColor_c)

%% 输出邻域生长法的结果
core=strel('disk',7);
u_plane_regiongrowing(depthColor_c, depthColor_c,9,core);

%% 输出去高光的结果
imshow(depthColor_c)
processed=u_basic_process(depthColor_c,3,core);