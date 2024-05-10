function [out,barrier_point,barrier_pos]=u_plane_regiongrowing(img_color, img_depth,nodesize,core,t_color)

[height,width,~]=size(img_color);

node_height=floor(height/nodesize);
node_width=floor(width/nodesize);
node_mid_index=floor(node_width/2);


nodes=cell(node_height, node_width);% nodes
% depth=zeros(node_height, node_width);
edges=ones(node_height,2).*node_mid_index;
% node_last_row=cell(1, node_width);
% scan_lines=cell(node_height, node_width);
% a=zeros(node_height, node_width);
% c=zeros(node_height, node_width);

% t_color=0.90;% 像素相似度的阈值
% t_merge=0.99;% 直线相似度的阈值
% left_stop=0;
% right_stop=0;

% 划分为node
if rem(nodesize,2) % nodesize为奇数，取中间1行
    for i=1:node_height
        for j=1:node_width
            nodes{i,j}=img_color((i-1)*nodesize+1:min(i*nodesize,height), ...
                (j-1)*nodesize+1:min(j*nodesize,width),:);% 划分彩色图
%             %             depth(i,j)=mean(img_depth((i-1)*nodesize+1:min(i*nodesize,height), ...
%             %                                       (j-1)*nodesize+1:min(j*nodesize,width)),"all");
%             % 以均值形式存储深度
%             scan_lines{i,j}=[(j-1)*nodesize+3,j*nodesize-2;...
%                 (i-1)*nodesize+1+floor(nodesize/2),(i-1)*nodesize+1+floor(nodesize/2);...
%                 img_depth((i-1)*nodesize+1+floor(nodesize/2),(j-1)*nodesize+3),...
%                 img_depth((i-1)*nodesize+1+floor(nodesize/2),j*nodesize-2)];
%             % 拟合扫描线
%             %             fit_spatial_line = fit(scan_lines{i,j}(1:2,:)',scan_lines{i,j}(3,:)',spatial_line);
%             a(i,j)=(scan_lines{i,j}(3,2)-scan_lines{i,j}(3,1))/(scan_lines{i,j}(1,2)-scan_lines{i,j}(1,1));
%             c(i,j)=scan_lines{i,j}(3,2)-a(i,j)*scan_lines{i,j}(1,2);
        end
    end

else % nodesize为偶数，取中间2行
    for i=1:node_height
        for j=1:node_width
            nodes{i,j}=img_color((i-1)*nodesize+1:min(i*nodesize,height), ...
                (j-1)*nodesize+1:min(j*nodesize,width),:);
%             %             depth(i,j)=mean(img_depth((i-1)*nodesize+1:min(i*nodesize,height), ...
%             %                                       (j-1)*nodesize+1:min(j*nodesize,width)),"all");
%             scan_lines{i,j}=[(j-1)*nodesize+3,j*nodesize-2,(j-1)*nodesize+3,j*nodesize-2;... %x
%                 (i-1)*nodesize+nodesize/2,(i-1)*nodesize+nodesize/2,(i-1)*nodesize+nodesize/2+1,(i-1)*nodesize+nodesize/2+1;... %y
%                 img_depth((i-1)*nodesize+nodesize/2,(j-1)*nodesize+3),img_depth((i-1)*nodesize+nodesize/2,j*nodesize-2),...%z
%                 img_depth((i-1)*nodesize+nodesize/2+1,(j-1)*nodesize+3),img_depth((i-1)*nodesize+nodesize/2+1,j*nodesize-2)];
%             a1=(scan_lines{i,j}(3,2)-scan_lines{i,j}(3,1))/(scan_lines{i,j}(1,2)-scan_lines{i,j}(1,1));
%             a2=(scan_lines{i,j}(3,4)-scan_lines{i,j}(3,3))/(scan_lines{i,j}(1,4)-scan_lines{i,j}(1,3));
%             c1=scan_lines{i,j}(3,2)-a1*scan_lines{i,j}(1,2);
%             c2=scan_lines{i,j}(3,4)-a2*scan_lines{i,j}(1,4);
%             a(i,j)=(a1+a2)/2;
%             c(i,j)=(c1+c2)/2;
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
% left_stop=node_mid_index-edges(node_height,1);
% edges(node_height,2) = max(idx);
% right_stop=edges(node_height,2)-node_mid_index;

% 从最下面一行开始扫描
left_stop=half_line_left(nodes(node_height,1:node_mid_index),node_mid_index,t_color);
right_stop=half_line_right(nodes(node_height,node_mid_index:end),t_color);
edges(node_height,1)=node_mid_index-left_stop;
edges(node_height,2)=node_mid_index+right_stop;
% nodes{node_height,edges(node_height,1)}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
% nodes{node_height,edges(node_height,2)}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));

%历史版本，未调用函数
% node=nodes{node_height,node_mid_index};
% color_mid=reshape(mean(node,[1,2]),1,3);
% node_left=nodes{node_height,node_mid_index-1};
% node_right=nodes{node_height,node_mid_index+1};
% color_left=reshape(mean(node_left,[1,2]),1,3);
% color_right=reshape(mean(node_right,[1,2]),1,3);
% left=colorvalue(color_left,color_mid,t_color);
% right=colorvalue(color_right,color_mid,t_color);
% for j=2:floor(node_width/2)+1
%     if left&&~left_stop %左侧色彩相似且没有停止，向左生长
%         node_left=nodes{node_height,node_mid_index-j};
%         color_left=reshape(mean(node_left,[1,2]),1,3);
%         left=colorvalue(color_left,color_mid,t_color);
%         if node_mid_index-j==1% 到达边界，停止
%             left_stop=j;
%             edges(node_height,1)=node_mid_index-left_stop;
%         end
%     elseif ~left&&~left_stop %色彩不相似且没有停止，停止生长
%         left_stop=j;
%         edges(node_height,1)=node_mid_index-left_stop;
%     end
%     if right&&~right_stop %右侧色彩相似且没有停止，向右生长
%         node_right=nodes{node_height,node_mid_index+j};
%         color_right=reshape(mean(node_right,[1,2]),1,3);
%         right=colorvalue(color_right,color_mid,t_color);
%         if node_mid_index+j==node_width
%             right_stop=j;
%             edges(node_height,2)=node_mid_index+right_stop;
%         end
%     elseif ~right&&~right_stop
%         right_stop=j-1;
%         edges(node_height,2)=node_mid_index+right_stop;
%     end
% end
% 底部一行可视化
% nodes{node_height,edges(node_height,1)}=uint8(repmat(reshape([255,0,0],1,1,3),nodesize,nodesize));% 给左边生长出的节点上色
% nodes{node_height,edges(node_height,2)}=uint8(repmat(reshape([255,0,0],1,1,3),nodesize,nodesize));
% imshow(cell2mat(nodes));


% 然后从下往上进行遍历
% 取两个node为一个判别指标。例如左边，要求左边node不相似，右边node相似，才判断出边界。
for i=node_height-1:-1:1 % node_height-11
    new_mid_index=floor((edges(i+1,1)+edges(i+1,2))/2);% 以下方一行的中间作为中间
    node=nodes{i,new_mid_index};

    node_left=nodes(i,node_mid_index-left_stop:node_mid_index-left_stop+1);%分配邻域
    %     imshow(cell2mat(node_left));
    node_right=nodes(i,node_mid_index+right_stop-1:node_mid_index+right_stop);
    %     imshow(cell2mat(node_right));
    left=neighbor(node_left,node,t_color);
    right=neighbor(node_right,node,t_color);

    % 判断左侧
    if isequal(left,[1,1])%两个节点都在道路上，重新按行扫描
        node_left=nodes(i,1:new_mid_index);
        left_stop=half_line_left(node_left,new_mid_index,t_color);
        edges(i,1)=new_mid_index-left_stop;
        left_stop=node_mid_index-edges(i,1);
    elseif isequal(left,[1,0])||isequal(left,[0,0])
        % 两个都在标线上，或者右侧在标线上，向右生长
        while ~isequal(left,[0,1]) && left_stop>=-node_mid_index
            left_stop=left_stop-1;
            node_left=nodes(i,node_mid_index-left_stop:node_mid_index-left_stop+1);
            left=neighbor(node_left,node,t_color);
        end
        edges(i,1)=node_mid_index-left_stop;
    else% [0，1]左侧在标线上，右侧不在
        edges(i,1)=node_mid_index-left_stop;
    end


    % 判断右侧
    if isequal(right,[1,1])%两个节点都在道路上，重新按行扫描
        node_right=nodes(i,new_mid_index:end);
        right_stop=half_line_right(node_right,t_color);
        edges(i,2)=new_mid_index+right_stop;
        right_stop=edges(i,2)-node_mid_index;
    elseif isequal(right,[0,1])||isequal(right,[0,0])
        % 两个都在标线上，或者左边在右边不在，向左生长
        while ~isequal(right,[1,0]) && right_stop>=-node_mid_index
            right_stop=right_stop-1;
            node_right=nodes(i,node_mid_index+right_stop-1:node_mid_index+right_stop);
            right=neighbor(node_right,node,t_color);
        end
        edges(i,2)=node_mid_index+right_stop;
    else%[1,0]的情况,左边不在标线上，右边在
        edges(i,2)=node_mid_index+right_stop;
    end
%     % 一行的扫描结果可视化
%     nodes{i,edges(i,1)}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
%     nodes{i,edges(i,2)}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
%     imshow(cell2mat(nodes));

end
% figure(2);imshow(cell2mat(nodes));title("邻域生长法线特征检测结果");
y=(1:node_height)';% y坐标以node最上边算
xl=edges(:,1);% x坐标以node最左边算
xr=edges(:,2);

% for i=1:node_height
%     nodes{i,edges(i,1)}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
%     nodes{i,edges(i,2)}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
% end
% imshow(cell2mat(nodes));

% 根据空间特征重新定义边界
% t_line=0.89;
% for i=floor(node_height/2):-1:1 % 根据空间特性重新定义边界
%     mid_slope=a(i,floor((edges(i,1)+edges(i,2))/2));
%     mid_intercept=c(i,floor((edges(i,1)+edges(i,2))/2));
%     for j=floor((edges(i,1)+edges(i,2))/2):-1:edges(i,1)
%         if ~m_line(a(i,j),mid_slope,c(i,j),mid_intercept,t_line)
%             edges(i,1)=j;
%             break;
%         end
%     end
%     for j=floor((edges(i,1)+edges(i,2))/2):1:edges(i,2)
%         if ~m_line(a(i,j),mid_slope,c(i,j),mid_intercept,t_line)
%             edges(i,2)=j;
%             break;
%         end
%     end
% end

% edges(16,2)=45;
% edges(17,2)=45;
% edges(18,2)=43;
% edges(19,2)=43;
% edges(20,2)=43;
% edges(21,2)=44;
% edges(22,2)=45;
% for i=1:node_height
%     nodes{i,edges(i,1)}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
%     nodes{i,edges(i,2)}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
% end
% imshow(cell2mat(nodes));


% 避障算法
barrier=img_depth;
barrier(floor(nodesize*node_height/3):end,:)=0;
barrier(barrier(:)>800)=0;
% barrier(barrier(:)>0)=1;
% out=cast(barrier,"logical");
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
    if rightmost_point<width/2 %在很近的距离上对齐失效，手动平移
         barrier_point_x= barrier_point_x-3*nodesize;
%     elseif leftmost_point>width/2
%         barrier_point_x= barrier_point_x+3*nodesize;
    end
    barrier_point_y=repmat((node_height-1).*nodesize+floor(nodesize/2)+1,1,size(barrier_point_x,2));

    barrier_pos=[leftmost_point,rightmost_point];
    out=[(y-1).*nodesize+1+floor(nodesize/2),(xl-1).*nodesize+1+floor(nodesize/2),(xr-1).*nodesize+1+floor(nodesize/2)];
    out=out(((end-10):end),:);
    barrier_point=[barrier_point_y',barrier_point_x',barrier_point_x'];
end
% 将边缘可视化
% img=insertMarker(img_color,[out(:,2),out(:,1)],'o','Color','green');
% img=insertMarker(img,[out(:,3),out(:,1)],'o','Color','magenta');
% imshow(img);

% % 在深度图上检验对应的区域是否为平面
% indexl=(xl-1).*nodesize+y;
% indexr=(xr-1).*nodesize+y;
% zl=depth(indexl);
% zr=depth(indexr);
%
%

% % 尝试合并节点，如果相邻的空间直线参数类似，则可以合并对应的节点
% xl=edges(:,1);% x坐标以node最左边算
% xr=edges(:,2);
% for i=1:node_height
%     aa(i)=mean(a(i,edges(i,1):edges(i,2)));
%     bb(i)=mean(c(i,edges(i,1):edges(i,2)));
% end
% xy=[];
% firstmerge=1;
% for k=2:node_height-1
%     m1=[aa(k-1),-1];
%     m2=[aa(k),-1];
%     m3=[aa(k+1),-1];
%     n1=[bb(k+1)-bb(k),-1];
%     n2=[bb(k)-bb(k-1),-1];
%     n3=[(bb(k+1)-bb(k-1))./2,-1];
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
% xy=(xy-1).*nodesize+floor(nodesize/2)*1.5;
% resultImage=img_color;
% 
% for k=1:20%size(xy,1)
%     resultImage = insertShape(resultImage, 'Line', xy(k,:), 'Color', 'green', 'LineWidth', nodesize+1);
% end

% figure(2); imshow(resultImage);


end

%% 计算像素相似程度，包含色彩相似和亮度相似
function output=colorvalue(rgb1,rgb2,t_color)
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
out=-(colorsin+valuesin)/2+1;
if out>=t_color
    output=1;
else
    output=0;
end
end

%% 在邻域中生长
function out=neighbor(nodes,node_mid,t_color)

color_mid=reshape(mean(node_mid,[1,2]),1,3);
node1=nodes{1,1};
node2=nodes{1,2};
% node3=nodes{1,3};
% node4=nodes{1,4};

color1=reshape(mean(node1,[1,2]),1,3);
color2=reshape(mean(node2,[1,2]),1,3);
% color3=reshape(mean(node3,[1,2]),1,3);
% color4=reshape(mean(node4,[1,2]),1,3);

out=[colorvalue(color1,color_mid,t_color),colorvalue(color2,color_mid,t_color)];
%             colorvalue(color3,color_mid,t_color),colorvalue(color4,color_mid,t_color)];
end

%% 在左边一半生长
function out=half_line_left(nodes,node_mid_index,t_color)
node=nodes{1,node_mid_index};
color_mid=reshape(mean(node,[1,2]),1,3);
node_left=nodes{1,node_mid_index-1};
color_left=reshape(mean(node_left,[1,2]),1,3);
left=colorvalue(color_left,color_mid,t_color);
left_stop=0;
for j=2:node_mid_index-1
    if left %左侧色彩相似且没有停止，向左生长
        node_left=nodes{1,node_mid_index-j};
        color_left=reshape(mean(node_left,[1,2]),1,3);
        left=colorvalue(color_left,color_mid,t_color);
        if node_mid_index-j==1% 到达边界，停止
            left_stop=j;
            break;
        end
    elseif ~left %色彩不相似，停止生长
        left_stop=j-1;
        break;
    end
end
out=left_stop;
end

%% 在右边一半生长
function out=half_line_right(nodes,t_color)
node=nodes{1,1};
[~,w]=size(nodes);
color_mid=reshape(mean(node,[1,2]),1,3);
right_stop=0;
for j=2:w
    node_right=nodes{1,j};
    color_right=reshape(mean(node_right,[1,2]),1,3);
    right=colorvalue(color_right,color_mid,t_color);
    if right %右侧色彩相似且没有停止，向右生长
        if j==w% 到达边界，停止
            right_stop=j-1;
            break;
        end
    elseif ~right %色彩不相似，停止生长
        right_stop=j-1;
        break;
    end
end
out=right_stop;
end

%% 合并在一个平面上的行
function out=m_line(a1,a2,b1,b2,t)
if (a1==0) || (a2==0)
    if b1==0 || b2==0 %缺失信息
        out=1;
    else % 平行于成像平面
        a1=sqrt(2)/2;
        a2=sqrt(2)/2;
        b=normalize([b1,b2],'norm',2);
        b1=b(1);
        b2=b(2);
        mLine=(a1*a2+b1*b2)/(sqrt(a1^2+b1^2)*sqrt(a2^2+b2^2))*0.5+0.5;
        if mLine>=t
            out=1;
        else
            out=0;
        end
    end
else %一般情况
    a=normalize([a1,a2],'norm',2);
    b=normalize([b1,b2],'norm',2);
    a1=a(1);
    a2=a(2);
    b1=b(1);
    b2=b(2);
    mLine=(a1*a2+b1*b2)/(sqrt(a1^2+b1^2)*sqrt(a2^2+b2^2))*0.5+0.5;
    if mLine>=t
        out=1;
    else
        out=0;
    end
end
end