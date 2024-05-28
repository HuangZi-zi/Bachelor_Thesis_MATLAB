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
% color=imread("Resource\snapc30.png");%  % %
% color=imread("Resource\curl.jpg");
% color=imread("Resource\snapc2.png");
color=imread("Resource\road.jpg");
% color=imread("Resource\snapc25.png");
color=imresize(color,[375,667]);
%depthColor_c=fliplr(imcrop(color,[89 1 width_cd height_cd-1]));
depthColor_c=imcrop(color,[89 1 width_cd height_cd-1]);

% depth=imread("Resource\snapd25.png");
depth=imread("Resource\snapd2.png");
% figure(3);imshow(depth)
% depthColor_d=fliplr(imcrop(depth,[1 8 width_cd-1 height_cd-1]));
depthColor_d=imcrop(depth,[1 8 width_cd height_cd-1]);
figure(1);imshow(depthColor_c)
figure(2);imshow(depthColor_d)
% figure(3);imshow(depth)
%% 输出邻域生长法的结果
core=strel('disk',7);
[edges,barrier,barrier_pos]=u_plane_regiongrowing(depthColor_c, depthColor_d,9,core,0.95);

%% 输出去高光的结果
imshow(depthColor_c)
processed=u_basic_process(depthColor_c,3,core);

%% 输出从轨迹-差速的映射函数
x=-300:1:300;
y=200./(1+exp(-(1/30)*x))-100;
plot(x,y);xlabel("偏离中心的像素");ylabel("两轮差速");

%% 输出APF轨迹规划

% edges_and_barrier=union(edges,barrier,"rows");
[out,dir1,dir2,weigh]=u_APF(depthColor_c,edges,100);
imshow(out)
% % 画箭头图
% [h,w,~]=size(depthColor_c);
% nodesize=edges(2,1)-edges(1,1);
% n=2*size(edges,1);%障碍个数
% %初始化车的参数
% Xo=[floor(w/2),h];%起点位置
% k=2;%计算引力需要的增益系数
% m=400;%计算斥力的增益系数
% Po=100;%障碍影响距离，当障碍和车的距离大于这个距离时，斥力为0，即不受该障碍的影响
% 
% l=0.5;%步长
% J=1000;%循环迭代次数
% 
% %给出障碍和目标信息，对于图像信息，先x再y
% destination=(edges(4,2)+edges(4,3)+edges(5,2)+edges(5,3)+edges(6,2)+edges(6,3))/6;% 三行节点求平均
% Xsum=[destination,edges(5,1);edges(:,2), edges(:,1);edges(:,3), edges(:,1)];
% %这个向量是(n+1)*2维，其中第一个点是目标位置，剩下的都是障碍的位置。
% XXX=Xo;%XXX表示车的当前位置
% 
% %***************初始化结束，开始主体循环******************
% for i=1:h/nodesize %循环开始
%     XXX(2)=(i-1)*nodesize+nodesize/2;
%     for j=1:w/nodesize
%         XXX(1)=(j-1)*nodesize+nodesize/2;
%         disp(XXX(1));
%         X((i-1)*w/nodesize+j)=XXX(1);
%         Y((i-1)*w/nodesize+j)=XXX(2);
%         for ii=1:n+1   %计算物体和障碍物、目标点的向量
%             deltaX(ii)=Xsum(ii,1)-XXX(1);
%             deltaY(ii)=Xsum(ii,2)-XXX(2);
%             r(ii)=sqrt(deltaX(ii)^2+deltaY(ii)^2);
%         end
%         %目标点对路径点的引力
%         Fatx=k*deltaX(1);
%         Faty=k*deltaY(1);
%         %各个障碍物对路径点的斥力
%         for ii=1:n
%             if r(ii+1)>Po
%                 Frex(ii)=0;
%                 Frey(ii)=0;
%             else
%                 Frex(ii)=-m*deltaX(ii+1)/(r(ii+1)*r(ii+1));
%                 Frey(ii)=-m*deltaY(ii+1)/(r(ii+1)*r(ii+1));
%             end
%         end
%         %计算合力
%         Fsumx((i-1)*w/nodesize+j)=Fatx+sum(Frex);
%         Fsumy((i-1)*w/nodesize+j)=Faty+sum(Frey);
%     end
% end
% out=insertMarker(depthColor_c,Xsum(2:end,:),'s','Color','red');% 障碍点
% out=insertShape(out,"filled-circle",[Xsum(1,1),Xsum(1,2),5],"Color",[255,165,0],Opacity=1);% 目标点
% out=insertShape(out,"filled-circle",[Xo,8],"Color","blue");% 起点
% 
% imshow(out);hold on
% quiver(X(end-w/nodesize*11:end),Y(end-w/nodesize*11:end),Fsumx(end-w/nodesize*11:end),Fsumy(end-w/nodesize*11:end),'g',LineWidth=1);
% %quiver(X(1:4:end),Y(1:4:end),Fsumx(1:4:end),Fsumy(1:4:end),'g',LineWidth=1);

%% 绘制阶跃响应
open('Resource/corrected.fig');%为文件名
handle = findobj(gca,'Type','line');%获取曲线的handle，如果图中有多条曲线，handle为一个数组
xdata = get(handle,'XData');
ydata = get(handle,'YData');
%%xlswrite('xdata.xlsx',xdata); %%我试过了这两行代码，保存下来的是空的xlsx文件，因此未采用
%%xlswrite('ydata.xlsx',ydata);
% figure(1);hold on
% plot(xdata,ydata)
xx=cell2mat(xdata(1:156)');
xx=xx-xx(end);
yy=cell2mat(ydata(1:156)');
temp=yy(1:85);
temp(temp(:)<-0.1)=temp(temp(:)<-0.1)+0.1;
temp=medfilt1(temp,8);
yy(1:85)=temp;
figure(3);
%hold on
plot(xx-4.3,yy,'.');
%figure(3);plot(xx,yy,'.');xlabel("时间");ylabel("预期两轮差速");

%% 绘制Hough变换结果
edge_img=u_edge(depthColor_c);
laneLines=u_line_hough(depthColor_c, edge_img);

%% 对比一次二维中值和两次一维中值
img=imread("Resource\oldman1.jpg");
img=im2gray(img);
subplot 221;imshow(img);title("origin",'FontName','Times New Roman');
img=imnoise(img,'gaussian',0.002);
img=imnoise(img,'salt & pepper',0.02);
img=im2double(img);
h=480;
w=640;

subplot 222;imshow(img);title("added noise",'FontName','Times New Roman');
% 中值滤波
filter_size=5;
% hist_eq_double = im2double(hist_eq);
% x_dir=hist_eq_double(:)';
% tic;
x_dir=img(:)';
fil_x=medfilt1(x_dir,filter_size);
fil_x_re=reshape(fil_x,h,w)';
y_dir=fil_x_re(:)';
fil_xy=medfilt1(y_dir,filter_size);%2次1维中值0.03024628s
fil_xy_re=reshape(fil_xy,w,h)';
% toc;
% disp(['运行时间: ',num2str(toc)]);
tic
filted=medfilt2(img,[filter_size,filter_size]);%1次2维中值0.18264729s
toc
disp(['运行时间: ',num2str(toc)]);
subplot 223;imshow(fil_xy_re);title("one-dimension mean-fil twice",'FontName','Times New Roman');
subplot 224;imshow(filted);title("two-dimensions mean-fil once",'FontName','Times New Roman');

%% 对比sobel和canny算子
img=imread("Resource\oldman1.jpg");
img=im2gray(img);
% subplot 221;imshow(img);title("origin");
img=imnoise(img,'gaussian',0.002);
img=imnoise(img,'salt & pepper',0.02);
img=im2double(img);
h=480;
w=640;

subplot 221;imshow(img);title("added noise",'FontName','Times New Roman');
% 中值滤波
filter_size=5;
% hist_eq_double = im2double(hist_eq);
% x_dir=hist_eq_double(:)';
% tic;
x_dir=img(:)';
fil_x=medfilt1(x_dir,filter_size);
fil_x_re=reshape(fil_x,h,w)';
y_dir=fil_x_re(:)';
fil_xy=medfilt1(y_dir,filter_size);%2次1维中值0.03024628s
fil_xy_re=reshape(fil_xy,w,h)';
subplot 222;imshow(fil_xy_re);title("filtered",'FontName','Times New Roman');

edgesobel = edge(fil_xy_re, 'sobel');
edgecanny = edge(fil_xy_re, "canny");
subplot 223;imshow(edgesobel);title("edge-Sobel",'FontName','Times New Roman');
subplot 224;imshow(edgecanny);title("edge-Canny",'FontName','Times New Roman');

%% 曲线测试情况
img=imread("Resource\snapc25.png");
% edge_img=u_edge(img);
% u_line_hough(img,edge_img);
u_plane_regiongrowing(img,img);

%% 人工势场法
img=imread("Resource\snapc30.png");
core=strel('disk',7);
edges=u_plane_regiongrowing(img,img,9,core);
[out,dir_this_time,dir_next_time,weigh]=u_APF(img,edges);
imshow(out);
disp(dir);