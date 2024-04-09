function [out,barrier_point,barrier_pos]=u_plane_regiongrowing(img_color, img_depth,nodesize,core)

[height,width,~]=size(img_color);

node_height=floor(height/nodesize);
node_width=floor(width/nodesize);
node_mid=floor(node_width/2);


nodes=cell(node_height, node_width);% nodes
% depth=zeros(node_height, node_width);
edges=ones(node_height,2).*node_mid;
% node_last_row=cell(1, node_width);
scan_lines=cell(node_height, node_width);
a=zeros(node_height, node_width);
c=zeros(node_height, node_width);

t_color=0.90;% 像素相似度的阈值
% t_merge=0.8;% 直线相似度的阈值
left_stop=0;
right_stop=0;

% 划分为node
if rem(nodesize,2) % nodesize为奇数，取中间1行
    for i=1:node_height
        for j=1:node_width
            nodes{i,j}=img_color((i-1)*nodesize+1:min(i*nodesize,height), ...
                                 (j-1)*nodesize+1:min(j*nodesize,width),:);
%             depth(i,j)=mean(img_depth((i-1)*nodesize+1:min(i*nodesize,height), ...
%                                       (j-1)*nodesize+1:min(j*nodesize,width)),"all");
            scan_lines{i,j}=[(j-1)*nodesize+3,j*nodesize-2;...
                            (i-1)*nodesize+1+floor(nodesize/2),(i-1)*nodesize+1+floor(nodesize/2);...
                            img_depth((i-1)*nodesize+1+floor(nodesize/2),(j-1)*nodesize+3),...
                            img_depth((i-1)*nodesize+1+floor(nodesize/2),j*nodesize-2)];
%             fit_spatial_line = fit(scan_lines{i,j}(1:2,:)',scan_lines{i,j}(3,:)',spatial_line);
            a(i,j)=(scan_lines{i,j}(3,2)-scan_lines{i,j}(3,1))/(scan_lines{i,j}(1,2)-scan_lines{i,j}(1,1));
            c(i,j)=scan_lines{i,j}(3,2)-a(i,j)*scan_lines{i,j}(1,2);
        end
    end
    
else % nodesize为偶数，取中间2行
    for i=1:node_height
        for j=1:node_width
            nodes{i,j}=img_color((i-1)*nodesize+1:min(i*nodesize,height), ...
                                 (j-1)*nodesize+1:min(j*nodesize,width),:);
%             depth(i,j)=mean(img_depth((i-1)*nodesize+1:min(i*nodesize,height), ...
%                                       (j-1)*nodesize+1:min(j*nodesize,width)),"all");
            scan_lines{i,j}=[(j-1)*nodesize+3,j*nodesize-2,(j-1)*nodesize+3,j*nodesize-2;... %x
                             (i-1)*nodesize+nodesize/2,(i-1)*nodesize+nodesize/2,(i-1)*nodesize+nodesize/2+1,(i-1)*nodesize+nodesize/2+1;... %y
                             img_depth((i-1)*nodesize+nodesize/2,(j-1)*nodesize+3),img_depth((i-1)*nodesize+nodesize/2,j*nodesize-2),...%z
                             img_depth((i-1)*nodesize+nodesize/2+1,(j-1)*nodesize+3),img_depth((i-1)*nodesize+nodesize/2+1,j*nodesize-2)];
            a1=(scan_lines{i,j}(3,2)-scan_lines{i,j}(3,1))/(scan_lines{i,j}(1,2)-scan_lines{i,j}(1,1));
            a2=(scan_lines{i,j}(3,4)-scan_lines{i,j}(3,3))/(scan_lines{i,j}(1,4)-scan_lines{i,j}(1,3));
            c1=scan_lines{i,j}(3,2)-a1*scan_lines{i,j}(1,2);
            c2=scan_lines{i,j}(3,4)-a2*scan_lines{i,j}(1,4);
            a(i,j)=(a1+a2)/2;
            c(i,j)=(c1+c2)/2;
        end
    end
end

% 以节点颜色为类，对最底部行进行聚类分析
% node_last_row=nodes(node_height,:);% 取出最下面一行节点
% image_array = cat(4, node_last_row{:});% 4维矩阵。装载这一行节点的三个通道。[node_size,node_size,channel,node_width]
% mean_colors = squeeze(mean(mean(image_array, 1), 2));%装载了一行节点三个通道分别的平均值[channel, node_width]
% 
% Z=linkage(mean_colors','average','chebychev');
% %figure();dendrogram(Z);
% clusters=cluster(Z,'cutoff',50,'criterion','distance')';
% clusters(clusters(:)~=1)=0;
% idx = find(clusters);% 找出c中非零元素的index
% 
% edges(node_height,1) = min(idx);
% left_stop=node_mid-edges(node_height,1);
% edges(node_height,2) = max(idx);
% right_stop=edges(node_height,2)-node_mid;

% 历史版本：从最下面中间开始遍历，再逐行向上
node=nodes{node_height,node_mid};
color_mid=reshape(mean(node,[1,2]),1,3);
node_left=nodes{node_height,node_mid-1};
node_right=nodes{node_height,node_mid+1};
color_left=reshape(mean(node_left,[1,2]),1,3);
color_right=reshape(mean(node_right,[1,2]),1,3);
left=colorvalue(color_left,color_mid);
right=colorvalue(color_right,color_mid);
%     nodes{node_height,node_mid}=uint8(repmat(reshape([255,0,0],1,1,3),nodesize,nodesize));% 给中间node上色
for j=2:floor(node_width/2)+1
    if left>t_color&&~left_stop %左侧色彩相似且没有停止，向左生长
%         nodes{node_height,node_mid-j}=uint8(repmat(reshape([255,0,0],1,1,3),nodesize,nodesize));% 给左边生长出的节点上色
        node_left=nodes{node_height,node_mid-j};
        color_left=reshape(mean(node_left,[1,2]),1,3);
        left=colorvalue(color_left,color_mid);
        if node_mid-j==1% 到达边界，停止
            left_stop=j;
            edges(node_height,1)=node_mid-left_stop;
        end
    elseif left<=t_color&&~left_stop %色彩不相似且没有停止，停止生长
        left_stop=j;
        edges(node_height,1)=node_mid-left_stop;
    end
    if right>t_color&&~right_stop %右侧色彩相似且没有停止，向右生长
%         nodes{node_height,node_mid+j}=uint8(repmat(reshape([255,0,0],1,1,3),nodesize,nodesize));
        node_right=nodes{node_height,node_mid+j};
        color_right=reshape(mean(node_right,[1,2]),1,3); 
        right=colorvalue(color_right,color_mid);
        if node_mid+j==node_width
            right_stop=j;
            edges(node_height,2)=node_mid+right_stop;
        end
    elseif right<=t_color&&~right_stop
        right_stop=j;
        edges(node_height,2)=node_mid+right_stop;
    end
end
% 底部一行可视化
% imshow(cell2mat(nodes));



% 然后从下往上进行遍历
for i=node_height-1:-1:floor(node_height/2)
    new_mid=floor((edges(i+1,1)+edges(i+1,2))/2);% 以下方一行的中间作为中间
    node=nodes{i,new_mid};
%     node=nodes{node_height-i,node_mid};
    node_left=nodes{i,node_mid-left_stop};% 向上生长一格
    node_right=nodes{i,node_mid+right_stop};
%     nodes{node_m-i,node_mid+right_stop-1}=uint8(repmat(reshape([0,255,0],1,1,3),10,10));
%     nodes{node_m-i,node_mid-left_stop+1}=uint8(repmat(reshape([0,0,255],1,1,3),10,10));
%     imshow(cell2mat(nodes));
    color_mid=reshape(mean(node,[1,2]),1,3);
    color_left=reshape(mean(node_left,[1,2]),1,3);
    color_right=reshape(mean(node_right,[1,2]),1,3);
    left=colorvalue(color_left,color_mid);
    right=colorvalue(color_right,color_mid);

    if left>t_color% 最左侧与中心相似，向左侧生长
        while left>t_color && left_stop<floor(node_width/2-1)
            left_stop=left_stop+1;
            node_left=nodes{i,node_mid-left_stop};
            color_left=reshape(mean(node_left,[1,2]),1,3);
            left=colorvalue(color_left,color_mid);
        end
%         nodes{i,node_mid-left_stop+1}=uint8(repmat(reshape([0,0,255],1,1,3),nodesize,nodesize));
        edges(i,1)=node_mid-left_stop+1;
    else% 最左侧与中心不相似，向右侧生长
        while left<t_color && left_stop>0
            left_stop=left_stop-1;
            node_left=nodes{i,node_mid-left_stop};
            color_left=reshape(mean(node_left,[1,2]),1,3);
            left=colorvalue(color_left,color_mid);
        end
%         nodes{i,node_mid-left_stop-1}=uint8(repmat(reshape([0,0,255],1,1,3),nodesize,nodesize));
        edges(i,1)=node_mid-left_stop-1;
    end

    if right>t_color% 最右侧与中心相似，向右侧生长
        while right>t_color && right_stop<floor(node_width/2)
            right_stop=right_stop+1;
            node_right=nodes{i,node_mid+right_stop};
            color_right=reshape(mean(node_right,[1,2]),1,3);
            right=colorvalue(color_right,color_mid);
        end
%         nodes{i,node_mid+right_stop-1}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
        edges(i,2)=node_mid+right_stop-1;
    else% 最右侧与中心不相似，向左侧扩展
        while right<t_color && right_stop>0
            right_stop=right_stop-1;
            node_right=nodes{i,node_mid+right_stop};
            color_right=reshape(mean(node_right,[1,2]),1,3);
            right=colorvalue(color_right,color_mid);
        end
%         nodes{i,node_mid+right_stop+1}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
        edges(i,2)=node_mid+right_stop+1;
    end
%      imshow(cell2mat(nodes));
end
% figure(2);imshow(cell2mat(nodes));title("邻域生长法线特征检测结果");
y=(1:node_height)';% y坐标以node最上边算
xl=edges(:,1);% x坐标以node最左边算
xr=edges(:,2);

% 根据空间特征重新定义边界
t_a=5;
t_c=100;
for i=node_height:-1:floor(node_height/2) % 根据空间特性重新定义边界
    mid_slope=a(i,floor((edges(i,1)+edges(i,2))/2));
    mid_intercept=c(i,floor((edges(i,1)+edges(i,2))/2));
    for j=floor((edges(i,1)+edges(i,2))/2):-1:edges(i,1)
        if abs(a(i,j)-mid_slope)>t_a || abs(c(i,j)-mid_intercept)>t_c
            edges(i,1)=j;
            break;
        end
    end
    for j=floor((edges(i,1)+edges(i,2))/2):1:edges(i,2)
        if abs(a(i,j)-mid_slope)>t_a || abs(c(i,j)-mid_intercept)>t_c
            edges(i,2)=j;
            break;
        end
    end
end

% 避障算法
barrier=img_depth;
barrier(floor(nodesize*node_height/3):end,:)=0;
barrier(barrier(:)>800)=0;
barrier=imbinarize(barrier);

barrier=imerode(barrier,core);
% figure(3);imshow(barrier);
[~, cols] = find(barrier);

if isempty(cols) % 没有障碍物
    barrier_pos=[width,1];
    out=[(y-1).*nodesize+1+floor(nodesize/2),(xl-1).*nodesize+1+floor(nodesize/2),(xr-1).*nodesize+1+floor(nodesize/2)];
    out=out(((end-10):end),:);
    barrier_point=[];
else % 有障碍物
    leftmost_point = min(cols);
    rightmost_point = max(cols);
    barrier_point_x=leftmost_point:nodesize:rightmost_point;
    barrier_point_y=repmat((node_height-3).*nodesize+6,1,size(barrier_point_x,2));

    barrier_pos=[leftmost_point,rightmost_point];
    out=[(y-1).*nodesize+1+floor(nodesize/2),(xl-1).*nodesize+1+floor(nodesize/2),(xr-1).*nodesize+1+floor(nodesize/2)];
    out=out(((end-10):end),:);
    barrier_point=[barrier_point_y',barrier_point_x',barrier_point_x'];
end



% % 在深度图上检验对应的区域是否为平面
% indexl=(xl-1).*nodesize+y;
% indexr=(xr-1).*nodesize+y;
% zl=depth(indexl);
% zr=depth(indexr);
% 
% 
% % 在深度图上拟合，一行节点对应一条空间直线
% a=(zl-zr)./(xl-xr);
% b=zl-a.*xl;
% 
% % 尝试合并节点，如果相邻的空间直线参数类似，则可以合并对应的节点
% xy=[];
% firstmerge=1;
% for k=2:node_m-1
%     m1=[a(k-1),-1];
%     m2=[a(k),-1];
%     m3=[a(k+1),-1];
%     n1=[b(k+1)-b(k),-1];
%     n2=[b(k)-b(k-1),-1];
%     n3=[(b(k+1)-b(k-1))./2,-1];
%     sum_m=dot(m1,m2)/(norm(m1)*norm(m2))+dot(m1,m3)/(norm(m1)*norm(m3))+dot(m2,m3)/(norm(m2)*norm(m3));
%     sum_n=dot(n1,n2)/(norm(n1)*norm(n2))+dot(n1,n3)/(norm(n1)*norm(n3))+dot(n2,n3)/(norm(n2)*norm(n3));
%     measure(k)=(sum_m+sum_n)/6+0.5;
%     if firstmerge && measure(k)>t_merge 
%         merged(k)=1;
%         xy=[xl(k-1),y(k-1),xr(k-1),y(k-1);xl(k),y(k),xr(k),y(k);xl(k+1),y(k+1),xr(k+1),y(k+1)];
%         firstmerge=0;
%     elseif measure(k)>t_merge && ~firstmerge
%         merged(k)=1;
%         xy=union(xy,[xl(k-1),y(k-1),xr(k-1),y(k-1);xl(k),y(k),xr(k),y(k);xl(k+1),y(k+1),xr(k+1),y(k+1)],"rows");
%     end
% end
% xy=(xy-1).*nodesize+1;
% resultImage=img_color;
% for k=1:size(xy,1)
%     resultImage = insertShape(resultImage, 'Line', xy(k,:), 'Color', 'blue', 'LineWidth', nodesize);
% end

% figure(2); imshow(resultImage);

% 检查深度图是否存在障碍，如果存在，则将障碍物添加进边界信息中。

end

%% 计算像素相似程度，包含色彩相似和亮度相似
function output=colorvalue(rgb1,rgb2)
    rgb1 = im2double(rgb1(:));
    rgb2 = im2double(rgb2(:));
    v1=0.3*rgb1(1)+0.59*rgb1(2)+0.11*rgb1(3);
    v2=0.3*rgb2(1)+0.59*rgb2(2)+0.11*rgb2(3);
    
    N1 = norm(rgb1);
    N2 = norm(rgb2);
    if isequal(rgb1,rgb2)
        colorsin = 0;
    else
        colorsin = norm(cross(rgb1,rgb2))/(N1*N2);
    end
    if isequal(v1,v2)
        valuesin = 0;
    else
        valuesin= norm(cross([v1;128;0],[v2;128;0])) / (norm([v1;128;0]) * (norm([v2;128;0])));
    end
    output=-(colorsin+valuesin)/2+1;
end

%% 检查平面相似度 optimality measure
% function output=optimality()
% 
% end