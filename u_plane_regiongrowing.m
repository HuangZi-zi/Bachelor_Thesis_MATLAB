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

t_color=3;% 颜色偏差的阈值
t_merge=0.8;% 直线相似度的阈值
left_stop=0;
right_stop=0;

% 划分为node
for i=1:node_m
    for j=1:node_n
        nodes{i,j}=img_color((i-1)*nodesize+1:min(i*nodesize,M),(j-1)*nodesize+1:min(j*nodesize,N),:);
        depth(i,j)=mean(img_depth((i-1)*nodesize+1:min(i*nodesize,M),(j-1)*nodesize+1:min(j*nodesize,N)),"all");
    end
end

% 先遍历最下面一行
node=nodes{node_m,node_mid};
color_mid=reshape(mean(node,[1,2]),1,3);
for j=1:floor(node_n/2-1)
    node_left=nodes{node_m,node_mid-j};
    node_right=nodes{node_m,node_mid+j};
    color_left=reshape(mean(node_left,[1,2]),1,3);
    color_right=reshape(mean(node_right,[1,2]),1,3);
    left=colorangle(color_left,color_mid);
    right=colorangle(color_right,color_mid);
    if left<t_color
        nodes{node_m,node_mid-j}=uint8(repmat(reshape([255,0,0],1,1,3),nodesize,nodesize));
        if j==floor(node_n/2-1)
            edges(node_m,1)=node_mid-j;
            left_stop=j;
        end
    elseif ~left_stop
        edges(node_m,1)=node_mid-j;
        left_stop=j;
    end
    if right<t_color&&~right_stop
        nodes{node_m,node_mid+j}=uint8(repmat(reshape([255,0,0],1,1,3),nodesize,nodesize));
        if j==floor(node_n/2-1)
            edges(node_m,2)=node_mid+j;
            right_stop=j;
        end
    elseif ~right_stop
        edges(node_m,2)=node_mid+j;
        right_stop=j;
    end
end

% imshow(cell2mat(nodes));
% 然后从下往上进行遍历

for i=1:node_m-1
    node=nodes{node_m-i,node_mid};
    node_left=nodes{node_m-i,node_mid-left_stop+1};
    node_right=nodes{node_m-i,node_mid+right_stop-1};
    %nodes{node_m-i,node_mid+right_stop-1}=uint8(repmat(reshape([0,255,0],1,1,3),10,10));
    %imshow(cell2mat(nodes));
    color_mid=reshape(mean(node,[1,2]),1,3);
    color_left=reshape(mean(node_left,[1,2]),1,3);
    color_right=reshape(mean(node_right,[1,2]),1,3);
    left=colorangle(color_left,color_mid);
    right=colorangle(color_right,color_mid);

    if left<t_color% 最左侧与中心相似，向左侧扩展
        while left<t_color && left_stop<node_n/2-1
            left_stop=left_stop+1;
            node_left=nodes{node_m-i,node_mid-left_stop+1};
            color_left=reshape(mean(node_left,[1,2]),1,3);
            left=colorangle(color_left,color_mid);
        end
        nodes{node_m-i,node_mid-left_stop+1}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
        edges(node_m-i,1)=node_mid-left_stop+1;
    else% 最左侧与中心不相似，向右侧扩展
        while left>t_color && left_stop>0
            left_stop=left_stop-1;
            node_left=nodes{node_m-i,node_mid-left_stop+1};
            color_left=reshape(mean(node_left,[1,2]),1,3);
            left=colorangle(color_left,color_mid);
        end
        nodes{node_m-i,node_mid-left_stop+1}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
        edges(node_m-i,1)=node_mid-left_stop+1;
    end

    if right<t_color% 最右侧与中心相似，向右侧扩展
        while right<t_color && right_stop<node_n/2
            right_stop=right_stop+1;
            node_right=nodes{node_m-i,node_mid+right_stop-1};
            color_right=reshape(mean(node_right,[1,2]),1,3);
            right=colorangle(color_right,color_mid);
        end
        nodes{node_m-i,node_mid+right_stop-1}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
        edges(node_m-i,2)=node_mid+right_stop-1;
    else% 最右侧与中心不相似，向左侧扩展
        while right>t_color && right_stop>0
            right_stop=right_stop-1;
            node_right=nodes{node_m-i,node_mid+right_stop-1};
            color_right=reshape(mean(node_right,[1,2]),1,3);
            right=colorangle(color_right,color_mid);
        end
        nodes{node_m-i,node_mid+right_stop-1}=uint8(repmat(reshape([0,255,0],1,1,3),nodesize,nodesize));
        edges(node_m-i,2)=node_mid+right_stop-1;
    end
end

% imshow(cell2mat(nodes));

y=(1:node_m)';% y坐标以node最上边算
xl=edges(:,1);% x坐标以node最左边算
xr=edges(:,2);

% x=[xl;xr];
% y=repmat(y,2,1);
indexl=(xl-1).*nodesize+y;
indexr=(xr-1).*nodesize+y;
zl=depth(indexl);
zr=depth(indexr);

a=(zl-zr)./(xl-xr);
b=zl-a.*xl;

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
% figure(1); imshow(resultImage);
end
timeit