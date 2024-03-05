function u_plane_regiongrowing(img_color, img_depth)

img_color=imcrop(img_color,[0,0,512,370]);
% figure; imshow(img_color);
[M,N,channel]=size(img_color);
% 720,1080
nodesize=10;
node_m=floor(M/nodesize);
node_n=floor(N/nodesize);
node_mid=floor(node_n/2);

nodes=cell(node_m, node_n);% nodes
depth=zeros(node_m, node_n);
edges=ones(node_m,2).*node_mid;
node_last_row=cell(1, node_n);

t_color=0.95;% 像素相似度的阈值
t_merge=0.8;% 直线相似度的阈值
left_stop=0;
right_stop=0;

% 划分为node
for i=1:node_m-1
    for j=1:node_n
        nodes{i,j}=img_color((i-1)*nodesize+1:min(i*nodesize,M),(j-1)*nodesize+1:min(j*nodesize,N),:);
        depth(i,j)=mean(img_depth((i-1)*nodesize+1:min(i*nodesize,M),(j-1)*nodesize+1:min(j*nodesize,N)),"all");
    end
end
i=node_m;
for j=1:node_n
    nodes{i,j}=img_color((i-1)*nodesize+1:min(i*nodesize,M),(j-1)*nodesize+1:min(j*nodesize,N),:);
    depth(i,j)=mean(img_depth((i-1)*nodesize+1:min(i*nodesize,M),(j-1)*nodesize+1:min(j*nodesize,N)),"all");
    node_last_row{1,j}=nodes{i,j};
end


% 先遍历最下面一行
% 对最下面一行进行聚类分析
image_array = cat(4, node_last_row{:});
mean_colors_last_row = squeeze(mean(mean(image_array, 1), 2)); % Calculate mean along rows and columns, then remove singleton dimensions
Z=linkage(mean_colors_last_row','average','chebychev');
figure();dendrogram(Z);
c=cluster(Z,'cutoff',80,'criterion','distance');
% 以最下面一行为种子节点，依次向上进行邻域生长


% 历史版本：从最下面中间开始遍历，再逐行向上
% node=nodes{node_m,node_mid};
% color_mid=reshape(mean(node,[1,2]),1,3);
% for j=1:floor(node_n/2-1)
%     node_left=nodes{node_m,node_mid-j};
%     node_right=nodes{node_m,node_mid+j};
%     color_left=reshape(mean(node_left,[1,2]),1,3);
%     color_right=reshape(mean(node_right,[1,2]),1,3);
%     left=colorvalue(color_left,color_mid);
%     right=colorvalue(color_right,color_mid);
%     if left>t_color
%         nodes{node_m,node_mid-j}=uint8(repmat(reshape([255,0,0],1,1,3),nodesize,nodesize));
%         if j==floor(node_n/2-1)
%             edges(node_m,1)=node_mid-j;
%             left_stop=j;
%         end
%     elseif ~left_stop
%         edges(node_m,1)=node_mid-j;
%         left_stop=j;
%     end
%     if right>t_color&&~right_stop
%         nodes{node_m,node_mid+j}=uint8(repmat(reshape([255,0,0],1,1,3),nodesize,nodesize));
%         if j==floor(node_n/2-1)
%             edges(node_m,2)=node_mid+j;
%             right_stop=j;
%         end
%     elseif ~right_stop
%         edges(node_m,2)=node_mid+j;
%         right_stop=j;
%     end
% end
% 
% % imshow(cell2mat(nodes));
% % 然后从下往上进行遍历
% 
% for i=1:node_m-1
%     node=nodes{node_m-i,node_mid};
%     node_left=nodes{node_m-i,node_mid-left_stop+1};
%     node_right=nodes{node_m-i,node_mid+right_stop-1};
%     %nodes{node_m-i,node_mid+right_stop-1}=uint8(repmat(reshape([0,255,0],1,1,3),10,10));
%     %imshow(cell2mat(nodes));
%     color_mid=reshape(mean(node,[1,2]),1,3);
%     color_left=reshape(mean(node_left,[1,2]),1,3);
%     color_right=reshape(mean(node_right,[1,2]),1,3);
%     left=colorvalue(color_left,color_mid);
%     right=colorvalue(color_right,color_mid);
% 
%     if left>t_color% 最左侧与中心相似，向左侧扩展
%         while left>t_color && left_stop<node_n/2-1
%             left_stop=left_stop+1;
%             node_left=nodes{node_m-i,node_mid-left_stop+1};
%             color_left=reshape(mean(node_left,[1,2]),1,3);
%             left=colorvalue(color_left,color_mid);
%         end
%         nodes{node_m-i,node_mid-left_stop+1}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
%         edges(node_m-i,1)=node_mid-left_stop+1;
%     else% 最左侧与中心不相似，向右侧扩展
%         while left<t_color && left_stop>0
%             left_stop=left_stop-1;
%             node_left=nodes{node_m-i,node_mid-left_stop+1};
%             color_left=reshape(mean(node_left,[1,2]),1,3);
%             left=colorvalue(color_left,color_mid);
%         end
%         nodes{node_m-i,node_mid-left_stop+1}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
%         edges(node_m-i,1)=node_mid-left_stop+1;
%     end
% 
%     if right>t_color% 最右侧与中心相似，向右侧扩展
%         while right>t_color && right_stop<node_n/2
%             right_stop=right_stop+1;
%             node_right=nodes{node_m-i,node_mid+right_stop-1};
%             color_right=reshape(mean(node_right,[1,2]),1,3);
%             right=colorvalue(color_right,color_mid);
%         end
%         nodes{node_m-i,node_mid+right_stop-1}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
%         edges(node_m-i,2)=node_mid+right_stop-1;
%     else% 最右侧与中心不相似，向左侧扩展
%         while right<t_color && right_stop>0
%             right_stop=right_stop-1;
%             node_right=nodes{node_m-i,node_mid+right_stop-1};
%             color_right=reshape(mean(node_right,[1,2]),1,3);
%             right=colorvalue(color_right,color_mid);
%         end
%         nodes{node_m-i,node_mid+right_stop-1}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
%         edges(node_m-i,2)=node_mid+right_stop-1;
%     end
% end
% 
% % imshow(cell2mat(nodes));


y=(1:node_m)';% y坐标以node最上边算
xl=edges(:,1);% x坐标以node最左边算
xr=edges(:,2);

% 将深度图反投影到点云
% [U,V,cloud_array]=projectPointCloud(X_t, Y_t, d_img, fx_rgb, fy_rgb, cx_rgb, cy_rgb, t_stereo(3));

% x=[xl;xr];
% y=repmat(y,2,1);
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

figure(1); imshow(resultImage);
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
%% 将深度图反投影为点云
function [U, V, cloud_array] = projectPointCloud(X, Y, Z, fx_rgb, fy_rgb, cx_rgb, cy_rgb, z_min)
    [height, width] = size(X);

    % Project to image coordinates
    U = X ./ Z;
    V = Y ./ Z;
    U = U * fx_rgb + cx_rgb;
    V = V * fy_rgb + cy_rgb;

    % Reusing U as cloud index
    %U = V * width + U + 0.5;

    cloud_array = zeros(height * width, 3);

    it = 1;
    for r = 1:height
        for c = 1:width
            z = Z(r, c);
            u = U(r, c);
            v = V(r, c);
            if (z > z_min && u > 0 && v > 0 && u < width && v < height)
                id = floor(v) * width + u;
                cloud_array(id, 1) = X(r, c);
                cloud_array(id, 2) = Y(r, c);
                cloud_array(id, 3) = z;
            end
            it = it + 1;
        end
    end
end