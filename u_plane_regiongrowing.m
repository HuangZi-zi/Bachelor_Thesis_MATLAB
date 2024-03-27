function out=u_plane_regiongrowing(img_color, img_depth,nodesize)
% out=[y，xl，xr]
% img_color=imcrop(img_color,[0,0,510,370]);
% img_depth=imcrop(img_depth,[0,0,510,370]);
% figure; imshow(img_color);

[height,width,channel]=size(img_color);

M=height;
N=width;
node_m=floor(M/nodesize);
node_n=floor(N/nodesize);
node_mid=floor(node_n/2);

nodes=cell(node_m, node_n);% nodes
depth=zeros(node_m, node_n);
edges=ones(node_m,2).*node_mid;
node_last_row=cell(1, node_n);

t_color=0.90;% 像素相似度的阈值
t_merge=0.8;% 直线相似度的阈值
left_stop=0;
right_stop=0;

% 划分为node
if rem(nodesize,2) % nodesize为奇数，取中间1行
    for i=1:node_height
        for j=1:node_width
            nodes{i,j}=img_color((i-1)*nodesize+1:min(i*nodesize,height), ...
                                 (j-1)*nodesize+1:min(j*nodesize,width),:);
            depth(i,j)=mean(img_depth((i-1)*nodesize+1:min(i*nodesize,height), ...
                                      (j-1)*nodesize+1:min(j*nodesize,width)),"all");
            scan_lines{i,j}=[(j-1)*nodesize+1:j*nodesize;...
                            repmat((i-1)*nodesize+1+floor(nodesize/2),1,nodesize);...
                            img_depth((i-1)*nodesize+1+floor(nodesize/2),(j-1)*nodesize+1:j*nodesize)];
%             fit_spatial_line = fit(scan_lines{i,j}(1:2,:)',scan_lines{i,j}(3,:)',spatial_line);
        end
    end
else % nodesize为偶数，取中间2行
    for i=1:node_height
        for j=1:node_width
            nodes{i,j}=img_color((i-1)*nodesize+1:min(i*nodesize,height), ...
                                 (j-1)*nodesize+1:min(j*nodesize,width),:);
            depth(i,j)=mean(img_depth((i-1)*nodesize+1:min(i*nodesize,height), ...
                                      (j-1)*nodesize+1:min(j*nodesize,width)),"all");
            scan_lines{i,j}=[(j-1)*nodesize+1:j*nodesize,(j-1)*nodesize+1:j*nodesize;...
                            repmat((i-1)*nodesize+floor(nodesize/2),1,nodesize),repmat((i-1)*nodesize+1+floor(nodesize/2),1,nodesize);...
                            img_depth((i-1)*nodesize+floor(nodesize/2),(j-1)*nodesize+1:j*nodesize),img_depth((i-1)*nodesize+1+floor(nodesize/2),(j-1)*nodesize+1:j*nodesize)];
        end
    end
end

% 以节点颜色为类，对最底部行进行聚类分析
% node_last_row=nodes(node_m,:);% 取出最下面一行节点
% image_array = cat(4, node_last_row{:});% 4维矩阵。装载这一行节点的三个通道。[node_size,node_size,channel,node_n]
% mean_colors = squeeze(mean(mean(image_array, 1), 2));%装载了一行节点三个通道分别的平均值[channel, node_n]
% 
% % line_top=img_depth(i*nodesize-9,:);
% % line_bottom=img_depth(i*nodesize,:);
% % reshaped_topline=reshape(line_top,nodesize,node_n);
% % reshaoed_bottomline=reshape(line_bottom,nodesize,node_n);
% % spatial_line = fittype('a*x + b*y + c', 'coefficients', {'a', 'b', 'c'}, 'independent', {'x', 'y'}, 'dependent', 'z');
% 
% % for j = 1:51    
% %     x = (1:nodesize)+(j-1)*nodesize;
% %     y_top = repmat(i*nodesize-9,[1,10]);
% %     y_bottom = repmat(i*nodesize,[1,10]);
% %     z_top = reshaped_topline(:, j);
% % 
% %     f_top = fit([x, y_top], z_top, spatial_line);
% %     f_bottom = fit([x, y_bottom], z_bottom, spatial_line);
% %     line_coefficients(:, i) = [f.a; f.b; f.c];
% % end
% 
% Z=linkage(mean_colors','average','chebychev');
% %figure();dendrogram(Z);
% c=cluster(Z,'cutoff',50,'criterion','distance')';
% c(c(:)~=1)=0;
% idx = find(c);% 找出c中非零元素的index
% idx_diff = diff(idx);% 在index中取微分
% break_indices = find(idx_diff > 1);% 如果微分大于1，则说明不连续
% edges(node_m,1) = [idx(1), idx(break_indices+1)];
% left_stop=node_mid-edges(node_m,1);
% edges(node_m,2) = [idx(break_indices), idx(end)];
% right_stop=edges(node_m,2)-node_mid;

% 历史版本：从最下面中间开始遍历，再逐行向上
node=nodes{node_m,node_mid};
color_mid=reshape(mean(node,[1,2]),1,3);
for j=1:floor(node_n/2-1)
    node_left=nodes{node_m,node_mid-j};
    node_right=nodes{node_m,node_mid+j};
    color_left=reshape(mean(node_left,[1,2]),1,3);
    color_right=reshape(mean(node_right,[1,2]),1,3);
    left=colorvalue(color_left,color_mid);
    right=colorvalue(color_right,color_mid);
    nodes{node_m,node_mid}=uint8(repmat(reshape([255,0,0],1,1,3),nodesize,nodesize));% 给中间node上色
    if left>t_color&&~left_stop %色彩相似且没有停止，向左生长
        nodes{node_m,node_mid-j}=uint8(repmat(reshape([255,0,0],1,1,3),nodesize,nodesize));% 给左边生长出的节点上色
        if j==floor(node_n/2-1)% 到达边界，停止
            left_stop=j;
            edges(node_m,1)=node_mid-left_stop;
        end
    elseif left<=t_color&&~left_stop %色彩不相似且没有停止，停止生长
        left_stop=j;
        edges(node_m,1)=node_mid-left_stop;
    end
    if right>t_color&&~right_stop
        nodes{node_m,node_mid+j}=uint8(repmat(reshape([255,0,0],1,1,3),nodesize,nodesize));
        if j==floor(node_n/2-1)
            right_stop=j;
            edges(node_m,2)=node_mid+right_stop;
        end
    elseif right<=t_color&&~right_stop
        right_stop=j;
        edges(node_m,2)=node_mid+right_stop;
    end
end
% 底部一行可视化
% imshow(cell2mat(nodes));



% 然后从下往上进行遍历
for i=1:node_m-1
%     new_mid=floor((edges(node_m,1)+edges(node_m,2))/2);% 以下方一行的中间作为中间
%     node=nodes{node_m-i,new_mid};
    node=nodes{node_m-i,node_mid};
    node_left=nodes{node_m-i,node_mid-left_stop};% 向上生长一格
    node_right=nodes{node_m-i,node_mid+right_stop};
%     nodes{node_m-i,node_mid+right_stop-1}=uint8(repmat(reshape([0,255,0],1,1,3),10,10));
%     nodes{node_m-i,node_mid-left_stop+1}=uint8(repmat(reshape([0,0,255],1,1,3),10,10));
%     imshow(cell2mat(nodes));
    color_mid=reshape(mean(node,[1,2]),1,3);
    color_left=reshape(mean(node_left,[1,2]),1,3);
    color_right=reshape(mean(node_right,[1,2]),1,3);
    left=colorvalue(color_left,color_mid);
    right=colorvalue(color_right,color_mid);

    if left>t_color% 最左侧与中心相似，向左侧生长
        while left>t_color && left_stop<floor(node_n/2-1)
            left_stop=left_stop+1;
            node_left=nodes{node_m-i,node_mid-left_stop};
            color_left=reshape(mean(node_left,[1,2]),1,3);
            left=colorvalue(color_left,color_mid);
        end
        nodes{node_m-i,node_mid-left_stop+1}=uint8(repmat(reshape([0,0,255],1,1,3),nodesize,nodesize));
        edges(node_m-i,1)=node_mid-left_stop+1;
    else% 最左侧与中心不相似，向右侧生长
        while left<t_color && left_stop>0
            left_stop=left_stop-1;
            node_left=nodes{node_m-i,node_mid-left_stop};
            color_left=reshape(mean(node_left,[1,2]),1,3);
            left=colorvalue(color_left,color_mid);
        end
        nodes{node_m-i,node_mid-left_stop-1}=uint8(repmat(reshape([0,0,255],1,1,3),nodesize,nodesize));
        edges(node_m-i,1)=node_mid-left_stop-1;
    end

    if right>t_color% 最右侧与中心相似，向右侧生长
        while right>t_color && right_stop<floor(node_n/2)
            right_stop=right_stop+1;
            node_right=nodes{node_m-i,node_mid+right_stop};
            color_right=reshape(mean(node_right,[1,2]),1,3);
            right=colorvalue(color_right,color_mid);
        end
        nodes{node_m-i,node_mid+right_stop-1}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
        edges(node_m-i,2)=node_mid+right_stop-1;
    else% 最右侧与中心不相似，向左侧扩展
        while right<t_color && right_stop>0
            right_stop=right_stop-1;
            node_right=nodes{node_m-i,node_mid+right_stop};
            color_right=reshape(mean(node_right,[1,2]),1,3);
            right=colorvalue(color_right,color_mid);
        end
        nodes{node_m-i,node_mid+right_stop+1}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
        edges(node_m-i,2)=node_mid+right_stop+1;
    end
%     imshow(cell2mat(nodes));
end
% figure(2);imshow(cell2mat(nodes));title("邻域生长法线特征检测结果");
y=(1:node_height)';% y坐标以node最上边算
xl=edges(:,1);% x坐标以node最左边算
xr=edges(:,2);

% 避障算法
barrier=img_depth;
barrier(barrier(:)>10000)=0;
barrier=imbinarize(barrier);
barrier(floor(height/3):end,:)=0;
% imshow(barrier);
[~, cols] = find(barrier);
leftmost_point = min(cols);
rightmost_point = max(cols);
barrier_point_x=leftmost_point:nodesize:rightmost_point;
barrier_point_y=repmat((node_height-3).*nodesize+6,1,size(barrier_point_x,2));

out=[(y-1).*nodesize+1+floor(nodesize/2),(xl-1).*nodesize+1+floor(nodesize/2),(xr-1).*nodesize+1+floor(nodesize/2)];
out=out(((end-10):end),:);
out=union(out,[barrier_point_y',barrier_point_x',barrier_point_x'],"rows");


% 在深度图上检验对应的区域是否为平面
indexl=(xl-1).*nodesize+y;
indexr=(xr-1).*nodesize+y;
zl=depth(indexl);
zr=depth(indexr);


% 在深度图上拟合，一行节点对应一条空间直线
a=(zl-zr)./(xl-xr);
b=zl-a.*xl;

% 尝试合并节点，如果相邻的空间直线参数类似，则可以合并对应的节点
xy=[];
firstmerge=1;
for k=2:node_m-1
    m1=[a(k-1),-1];
    m2=[a(k),-1];
    m3=[a(k+1),-1];
    n1=[b(k+1)-b(k),-1];
    n2=[b(k)-b(k-1),-1];
    n3=[(b(k+1)-b(k-1))./2,-1];
    sum_m=dot(m1,m2)/(norm(m1)*norm(m2))+dot(m1,m3)/(norm(m1)*norm(m3))+dot(m2,m3)/(norm(m2)*norm(m3));
    sum_n=dot(n1,n2)/(norm(n1)*norm(n2))+dot(n1,n3)/(norm(n1)*norm(n3))+dot(n2,n3)/(norm(n2)*norm(n3));
    measure(k)=(sum_m+sum_n)/6+0.5;
    if firstmerge && measure(k)>t_merge 
        merged(k)=1;
        xy=[xl(k-1),y(k-1),xr(k-1),y(k-1);xl(k),y(k),xr(k),y(k);xl(k+1),y(k+1),xr(k+1),y(k+1)];
        firstmerge=0;
    elseif measure(k)>t_merge && ~firstmerge
        merged(k)=1;
        xy=union(xy,[xl(k-1),y(k-1),xr(k-1),y(k-1);xl(k),y(k),xr(k),y(k);xl(k+1),y(k+1),xr(k+1),y(k+1)],"rows");
    end
end
xy=(xy-1).*nodesize+1;
resultImage=img_color;
for k=1:size(xy,1)
    resultImage = insertShape(resultImage, 'Line', xy(k,:), 'Color', 'blue', 'LineWidth', nodesize);
end

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
function output=optimality()

end