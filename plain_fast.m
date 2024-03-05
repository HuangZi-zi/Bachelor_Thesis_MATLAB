clear;
img=imread("C:\Users\YawnFun\Downloads\rgbd_dataset_freiburg1_floor\depth\1305033539.743179.png");
rgb_img=imread("C:\Users\YawnFun\Downloads\rgbd_dataset_freiburg1_floor\rgb\1305033539.743168.png");
[M,N,channel]=size(img);
height=M;
width=N;
H=16;
W=16;

% img=double(im2gray(img));
img=cast(img,'double');
% img=histeq(img,256);
% img=img.*255;

% filter_size=9;
% x_dir=img(:)';
% fil_x=medfilt1(x_dir,filter_size);
% fil_x_re=reshape(fil_x,M,N)';
% y_dir=fil_x_re(:)';
% fil_xy=medfilt1(y_dir,filter_size);
% img=reshape(fil_xy,N,M)';

% figure();imshow(img);

Tmse=400;

% data=zeros(M*N,3);
% for i=1:M
%     for j=1:N
%         data((i-1)*N+j,:)=[i;j;img(i,j)];
%     end
% end

numNodesRows = floor(M / H);
numNodesCols = floor(N / W);

V_pixel=cell(numNodesRows, numNodesCols);% nodes
V_normal=zeros(numNodesRows, numNodesCols, 3);
V_angle=zeros(numNodesRows, numNodesCols);
V_plane=cell(numNodesRows, numNodesCols);
V_MSE=zeros(numNodesRows, numNodesCols);

V_merge=zeros(M,N,3);

merged=cell(1,numNodesRows*numNodesCols/2);

x_mean=zeros(numNodesRows, numNodesCols);
y_mean=zeros(numNodesRows, numNodesCols);
z_mean=zeros(numNodesRows, numNodesCols);

% 将深度图反投影为点云
% 校准数据
K_rgb=[5.4886723733696215e+02 0. 3.1649655835885483e+02;
    0. 5.4958402532237187e+02 2.2923873484682150e+02;
    0. 0. 1.];
dist_coeffs_rgb=[7.5390330809303213e-02 -2.8252164126047641e-02 ...
    -9.5892915347463477e-03 -4.6479123130854875e-04 ...
    -1.4030650062583919e-01];
K_ir=[5.7592685448804468e+02 0. 3.1515026356388171e+02;
    0. 5.7640791601093247e+02 2.3058580662101753e+02;
    0. 0. 1.];
dist_coeffs_ir=[-4.5096285298902479e-02 2.3953450975051985e-01 ...
    -9.5873062860576726e-03 -2.9427017385906013e-03 ...
    -1.0147872415919272e+00];
R_stereo=[9.9993556622059065e-01 -1.1131100247700660e-02 2.2275579414240366e-03;
    1.1108891478543052e-02 9.9989080121799856e-01 9.7456744121138719e-03;
    -2.3357947736726680e-03 -9.7203007620451660e-03 9.9995002865936788e-01];
t_stereo=[1.1497548441022023e+01; 3.5139088879273231e+01; 2.1887459420807019e+01];
fx_ir = K_ir(1,1);
fy_ir = K_ir(2,2);
cx_ir = K_ir(1,3);
cy_ir = K_ir(2,3);
fx_rgb = K_rgb(1,1);
fy_rgb = K_rgb(2,2);
cx_rgb = K_rgb(1,3);
cy_rgb = K_rgb(2,3);
% 预计算
X_pre=cast(zeros(height,width),'double');
Y_pre=cast(zeros(height,width),'double');
for r=0:1:height-1
    for c=0:1:width-1
        X_pre(r+1,c+1) = (c-cx_ir)/fx_ir;
        Y_pre(r+1,c+1) = (r-cy_ir)/fy_ir;
    end
end
X = X_pre.*img;
Y = Y_pre.*img;
X_t = (R_stereo(1,1))*X+(R_stereo(1,2))*Y+(R_stereo(1,3))*img + t_stereo(1);
Y_t = (R_stereo(2,1))*X+(R_stereo(2,2))*Y+(R_stereo(2,3))*img + t_stereo(2);
d_img = (R_stereo(3,1))*X+(R_stereo(3,2))*Y+(R_stereo(3,3))*img + t_stereo(3);
[U,V,cloud_array]=projectPointCloud(X_t, Y_t, d_img, fx_rgb, fy_rgb, cx_rgb, cy_rgb, t_stereo(3));
reshaped_point_cloud = reshape(cloud_array, [640, 480, 3]);
reshaped_point_cloud=pagetranspose(reshaped_point_cloud);

% v=struct('pixels', [], 'normal', [], 'plane', [],'v_MSE',[]);
E=zeros(H*W,4);% edge

% 主成分分析
for i=1:numNodesRows
    for j=1:numNodesCols
        V_pixel{i,j}=reshaped_point_cloud((i-1)*H+1:min(i*H,M),(j-1)*W+1:min(j*W,N),:);
%         % 每个node中的相对座标
%         points_r(:,1)=repelem(1:H,W)';%x
%         points_r(:,2)=repmat(1:W,1,W)';%y
%         % 归化为均值为0的序列后每个node中的座标
%         points_c(:,1)=points_r(:,1)-(1+W)/2;
%         points_c(:,2)=points_r(:,2)-(1+H)/2;
        points(:,1)=reshape(V_pixel{i,j}(:,:,1),[],1);%x
        points(:,2)=reshape(V_pixel{i,j}(:,:,2),[],1);%y
        points(:,3)=reshape(V_pixel{i,j}(:,:,3),[],1);%z

        x_mean(i,j)=mean(points(:,1));
        y_mean(i,j)=mean(points(:,2));
        z_mean(i,j)=mean(points(:,3));

        points_c(:,1)=points(:,1)-x_mean(i,j);
        points_c(:,2)=points(:,2)-y_mean(i,j);
        points_c(:,3)=points(:,3)-z_mean(i,j);

        %[V_normal{i,j}, V_plane{i,j}] = fitPlane_SVD(V_pixel{i,j});

%         var_x=std(x)^2;
%         var_y=std(y)^2;
%         var_z=std(z)^2;
%         cov_xy=cov(x,y);
%         cov_xz=cov(x,z);
%         cov_yz=cov(y,z);
%         sigma=[var_x,cov_xy,cov_xz;cov_xy,var_y,cov_yz;cov_xz,cov_yz,var_z];

        sigma=points_c'*points_c./(W*H);

        [eigenvector,eigenvalue]=eig(sigma);
        eigenvalue=diag(eigenvalue);
        [sorted_enigenvalue, index]=sort(eigenvalue);
        normal_tempt=eigenvector(index(1),:)./norm(eigenvector(index(1),:));
        V_normal(i,j,:)=reshape(normal_tempt, [1, 1, 3]);
        V_angle(i,j)=acos(dot(normal_tempt,[1,1,1])/(norm(normal_tempt)*norm([1,1,1])))/pi*180;
        %         gmag=imgradient(V_pixel{i,j});
        
        dx = gradient(points_c(:,1));
        dy = gradient(points_c(:,2));
        dz = gradient(points_c(:,3));
        gmag = sqrt(dx.^2+dy.^2+dz.^2);

        if  any([z_mean(i,j)<2,max(max(gmag))>200,sorted_enigenvalue(1)>Tmse]) %missing %discontinue %large mse
            V_MSE(i,j)=inf;
            V_normal(i,j)=inf;
            V_angle(i,j)=inf;
            V_pixel{i,j}=[];
            z_mean(i,j)=inf;
            E=union(E,[i,j-1,i,j ; i,j,i,j+1 ; i-1,j,i,j ; i,j,i+1,j],'rows');
        else
            V_MSE(i,j)=sorted_enigenvalue(1);
        end  
    end
end

% V_MSE_T=V_MSE';
% heap=(V_MSE_T(:))';
% [sorted_heap, index]=sort(heap);
% d=pdist(heap,@distfun);
% z=linkage(heap,'average','@distfun');

nodes(:,1)=repelem(1:numNodesRows,numNodesCols)';%x
nodes(:,2)=repmat(1:numNodesCols,1,numNodesRows)';%y
V_angle=V_angle';
nodes(:,3)=V_angle(:)';
z_mean=z_mean';
nodes(:,4)=z_mean(:)';

Z=linkage(nodes,'average','chebychev');
figure();dendrogram(Z);
c=cluster(Z,'cutoff',400,'criterion','distance');
% figure()
% gscatter(nodes(:,2),nodes(:,1),c);

% Define a colormap for visualization
colormap = rand(500, 3); % Adjust the size (10) based on the maximum number of clusters

% Loop through each node in the grid
for row = 1:numNodesRows
    for col = 1:numNodesCols
        % Get the cluster assignment for the current node
        currentNodeCluster = c((row - 1) * numNodesCols + col);
        
        % Get the color for the current cluster
        currentColor = colormap(currentNodeCluster, :);
        
        % Assign the color to the pixels of the current node in the output image
        V_merge((row - 1) * 16 + 1 : row * 16, (col - 1) * 16 + 1 : col * 16, :) = repmat(reshape(currentColor, [1, 1, 3]), [16,16,1]);
    end
end

% Display the resulting image
figure();imshow(V_merge);
figure();imshow(rgb_img);

% function distance = distfun(ZI,ZJ)
% [xrow,xcolumn]=size(ZI);
% [yrow,ycolumn]=size(ZJ);
% %初始化距离
% distance=zeros(yrow,1);
% if (xrow==1 && xcolumn==ycolumn)
%     for m=1:yrow
%         x1=ZI;%必须是行向量,不能是空向量
%         y1=ZJ(m,:);%必须是行向量,不能是空向量
%         b=( ~isnan(x1)) & (~isnan(y1)); %提取(x1,y1)中都不是nan的索引
%         A=[];
%         A(1,:)=x1(b);%必须是行向量,不能是空向量
%         A(2,:)=y1(b);%必须是行向量,不能是空向量
%         %计算距离
%         [~,distance(m,1)]=trymerge(A,100);
%     end
% end
% 
% end


% merge_cluster=1;
% 
% for i=1:size(heap,2)
%     current=index(1);
%     current_i=floor((current-1)/numNodesCols)+1;
%     current_j=current-(current_i-1)*numNodesCols;%rem(current,numNodesCols);
%     if current_i==1||current_i==numNodesRows||current_j==1||current_j==numNodesCols
%         continue
%     end
% 
%     if ~isempty(V_pixel{current_i,current_j})
%         points_current(:,1)=points_r(:,1)+H*(current_i-1);%x
%         points_current(:,2)=points_r(:,2)+W*(current_j-1);%y
%         points_current(:,3)=(V_pixel{current_i,current_j}(:));%z
% %         index(1)=[];
%         if isempty(merged{1,1})
%             merged{1,1}=points_current; 
%             V_merge(points_current(:,1),points_current(:,2))=merge_cluster;
%             index=setdiff(index,current,'stable');
%         else
%             for j=1:merge_cluster
%                 older=union(merged{1,j},points_current,"rows");
%                 [outcome,newMse]=trymerge(older,Tmse);
%                 if outcome
%                     merge_cluster=j;
%                     V_merge(points_current(:,1),points_current(:,2))=merge_cluster;
%                     merged{1,merge_cluster}=union(merged{1,merge_cluster},points_current,"rows");
%                     index=setdiff(index,current,'stable');
%                     break;
%                 elseif j==merge_cluster
%                     % 如果没有break，说明不能merge进之前的块里
%                     merge_cluster=merge_cluster+1;
%                     merged{1,merge_cluster}=points_current;
%                 end
%             end
% 
%         end
%     else
%         continue
%     end
%     
%     if ~isempty(V_pixel{current_i,current_j-1}) && ismember(current-1,index)
%         points_left(:,1)=points_r(:,1)+H*(current_i-1);%x
%         points_left(:,2)=points_r(:,2)+W*(current_j-1-1);%y
%         points_left(:,3)=(V_pixel{current_i,current_j-1}(:));%z
%         merged_left=union(points_current,points_left,"rows");
%         [outcome,newMse]=trymerge(merged_left,Tmse);
%         if outcome
%             V_merge(points_left(:,1),points_left(:,2))=merge_cluster;
%             merged{1,merge_cluster}=union(merged{1,merge_cluster},merged_left,"rows");
%             index=setdiff(index,current-1,'stable');
% %         else
% %             merge_cluster=merge_cluster+1;
% %             merged{1,merge_cluster}=points_current;
%         end
%     end
%     
%     if ~isempty(V_pixel{current_i,current_j+1}) && ismember(current+1,index)
%         points_right(:,1)=points_r(:,1)+H*(current_i-1);%x
%         points_right(:,2)=points_r(:,2)+W*(current_j-1+1);%y
%         points_right(:,3)=(V_pixel{current_i,current_j+1}(:));%z
%         merged_right=union(points_current,points_right,"rows");
%         [outcome,newMse]=trymerge(merged_right,Tmse);
%         if outcome
%             V_merge(points_right(:,1),points_right(:,2))=merge_cluster;
%             merged{1,merge_cluster}=union(merged{1,merge_cluster},merged_right,"rows");
%             index=setdiff(index,current+1,'stable');
% %         else
% %             merge_cluster=merge_cluster+1;
% %             merged{1,merge_cluster}=points_current;
%         end
%     end
% 
%     if ~isempty(V_pixel{current_i-1,current_j}) && ismember(current-numNodesCols,index)
%         points_up(:,1)=points_r(:,1)+H*(current_i-1-1);%x
%         points_up(:,2)=points_r(:,2)+W*(current_j-1);%y
%         points_up(:,3)=(V_pixel{current_i-1,current_j}(:));%z
%         merged_up=union(points_current,points_up,"rows");
%         [outcome,newMse]=trymerge(merged_up,Tmse);
%         if outcome
%             V_merge(points_up(:,1),points_up(:,2))=merge_cluster;
%             merged{1,merge_cluster}=union(merged{1,merge_cluster},merged_up,"rows");
%             index=setdiff(index,current-numNodesCols,'stable');
% %         else
% %             merge_cluster=merge_cluster+1;
% %             merged{1,merge_cluster}=points_current;
%         end
%     end
% 
%     if ~isempty(V_pixel{current_i+1,current_j}) && ismember(current+numNodesCols,index)
%         points_down(:,1)=points_r(:,1)+H*(current_i-1+1);%x
%         points_down(:,2)=points_r(:,2)+W*(current_j-1);%y
%         points_down(:,3)=(V_pixel{current_i+1,current_j}(:));%z
%         merged_down=union(points_current,points_down,"rows");
%         [outcome,newMse]=trymerge(merged_down,Tmse);
%         if outcome
%             V_merge(points_down(:,1),points_down(:,2))=merge_cluster;
%             merged{1,merge_cluster}=union(merged{1,merge_cluster},merged_down,"rows");
%             index=setdiff(index,current+numNodesCols,'stable');
% %         else
% %             merge_cluster=merge_cluster+1;
% %             merged{1,merge_cluster}=points_current;
%         end
%     end
%     imshow(uint8(V_merge.*100));
% %     right=current+1;
% %     left=current-1;
% %     up=current-numNodesCols;
% %     down=current+numNodesCols;  
% end

%imshow(uint8(V_merge.*10));
    
% function [outcome,newMSE]=trymerge(points,Tmse)
% W=16;H=16;
% x=points(:,1);
% y=points(:,2);
% z=points(:,3);
% 
% pc(:,1)=x-mean(x);
% pc(:,2)=y-mean(y);
% pc(:,3)=z-mean(z);
% 
% sigma=pc'*pc./(W*H);
%         [eigenvector,eigenvalue]=eig(sigma);
%         eigenvalue=diag(eigenvalue);
%         [sorted_enigenvalue, index]=sort(eigenvalue);
%         if sorted_enigenvalue(1)<Tmse
%             outcome=true;
%             newMSE=sorted_enigenvalue(1);
%         else
%             outcome=false;
%             newMSE=inf;
%         end
% 
% end



%%

% for i=2:M/H
%     for j=2:N/W
%         if rejectedge(V_pixel{i,j-1},V_pixel{i,j},V_pixel{i,j+1},V_normal{i,j-1},V_normal{i,j+1})
%             E=union(E,[V(i,j-1), V(i,j);V(i,j), V(i,j+1)]);
%         elseif rejectedge(V{i-1,j},V{i,j},V{i+1,j})
%             E=union(E,[V(i-1,j) V(i,j);V(i,j) V(i+1,j)]);
%         end
%     end
% end



%% 



% %%
% function [output1,output2] = initgraph(img)
% [M,N,channel]=size(img);
% H=16;
% W=16;
% 
% numNodesRows = floor(M / H);
% numNodesCols = floor(N / W);
% 
% V_pixel=cell(numNodesRows, numNodesCols);% nodes
% V_normal=cell(numNodesRows, numNodesCols);
% V_plane=cell(numNodesRows, numNodesCols);
% V_v_MSE=cell(numNodesRows, numNodesCols);
% % v=struct('pixels', [], 'normal', [], 'plane', [],'v_MSE',[]);
% E=[];% edge
% for i=1:numNodesRows
%     for j=1:numNodesCols
%         V_pixel{i,j}=img((i-1)*H+1:min(i*H,M),(j-1)*W+1:min(j*W,N));
%         [V_normal{i,j}, V_plane{i,j}, V_v_MSE{i,j}] = fitPlane_SVD(V_pixel{i,j});
%         if rejectnode(V_pixel{i,j},V_v_MSE{i,j})
%             V_pixel{i,j}=[];%zeros(H,W);
%             V_v_MSE{i,j}=inf;         
%         end
%     end
% end
% 
% for i=2:M/H
%     for j=2:N/W
%         if rejectedge(V_pixel{i,j-1},V_pixel{i,j},V_pixel{i,j+1},V_normal{i,j-1},V_normal{i,j+1})
%             E=union(E,[V(i,j-1), V(i,j);V(i,j), V(i,j+1)]);
%         elseif rejectedge(V{i-1,j},V{i,j},V{i+1,j})
%             E=union(E,[V(i-1,j) V(i,j);V(i,j) V(i+1,j)]);
%         end
%     end
% end
% output1=V;
% output2=E;
% end
% 
%%
% function output=rejectnode(input1,input2)
% % 1输入像素 2输入mse
% v=input1;
% [H,W,channel]=size(v);
% gmag=imgradient(v);
% Tmse=0.1;
% count_miss=0;
% for i=1:H
%     for j=1:W
%         if v(i,j)==0 % missing data
%             count_miss=count_miss+1;
%         elseif gmag(i,j)> 20 % depth-discontinuous
%             output=true;
%             return;
%         end
%     end
% end
% if count_miss>25
%     output=true;
%     return;
% elseif input2>Tmse
%     output=true;
%     return;
% else
%     output=false;
% end
% end

%%
% function output=rejectedge(va,vb,vc,na,nc)
% % 123输入像素 45输入法向量
% Tang=10;
% normal_va=na;
% normal_vc=nc;
% sigma = acos(dot(normal_va,normal_vc)/(norm(normal_va)*norm(normal_vc)));%弧度制
% normal_diff=abs(sigma/pi*180); %换算成角度
% 
% % if va==0||vb==0||vc==0
% if isempty(va)||isempty(vb)||isempty(vc)
%     output=true;
%     return;
% elseif normal_diff>Tang
%     output=true;
%     return;
% else
%     output=false;
% end
% end
% 
% %% 
% function [normal, planes] = fitPlane_SVD(v)
% % 功能：利用SVD拟合平面
% % 输入：data   - 原始数据(m*3)    
% % 输出：planes - 拟合所得平面参数 
% % v=input.pixels
% [w,h,channel]=size(v);
% data=zeros(w*h,3);
% for i=1:w
%     for j=1:h
%         data((i-1)*h+j,:)=[i;j;v(i,j)];
%     end
% end
% points = data(:,1:3);
% % 去质心
% M = points - repmat(mean(points),size(points,1),1);
% % 奇异值分解
% [~,~,V] = svd(M,0);
% % 最小特征值对应的特征向量为法向量
% normal = V(:,3)';
% % 平面参数标准化
% dtmp = mean(points*normal');
% planes(1:3) = normal'*sign(dtmp);
% planes(4) = -dtmp*sign(dtmp);
% 
% % x = 1:w;
% % y = 1:h;
% % [X,Y]= meshgrid (x,y);
% % ZFit = -(planes(4) + planes(1) * X + planes(2) * Y)/planes(3);
% % error=v-ZFit;
% % v_MSE=mse(error(:));
% end
% 
% % %%
% % % function outcome=MSE(input)
% % % v=input.pixels;
% % % if v==0
% % %     outcome=inf;
% % %     return;
% % % else 
% % %     [~,outcome]=fitPlane_SVD(input);
% % % end
% % % end
% % 
% % %%
% function ourput=AHCluster(V,E)
% Tmse=0.1;
% Tnum=160;
% Q=buildminMSEheap(V);% build a min-heap data structure
% B=[];
% PI=[];% 平面的方程
% while ~isempty(Q)
%     v=popmin(Q);
%     if v notin V
%         continue
%     end
%     u_best=[];
%     u_merge=[];
%     for i=1:size(E,2)
%         u_test=union(E(i,1),E(i,2));
%         if MSE(u_test)<MSE(u_merge)
%             u_best=E(i,1);
%             u_merge=u_test;
%         end
%     end
%     if MSE(u_merge)>=Tmse
%         if abs(v)>=Tnum
%             B=union(B,v);
%             norm=fitPlane_SVD(v);
%             dtmp = mean(points*norm');
%             planes(1:3) = norm'*sign(dtmp);
%             planes(4) = -dtmp*sign(dtmp);
%             PI=union(PI,planes);
%             %从E中删除E（v）
%             E=setdiff(E,E(v));
%             %从V中删除v
%             V=setdiff(V,v);
%         else
%             insert(Q,u_merge);
%             %E=union(E,)
%         end
%     end
% 
%     
% end
% end


%% 
% 新思路：用梯度来聚类平面
% clear;
% img=imread("Resource\depth.png");
% [M,N,channel]=size(img);
% img=im2double(img);
% img=histeq(img,256);
% 
% filter_size=9;
% x_dir=img(:)';
% fil_x=medfilt1(x_dir,filter_size);
% fil_x_re=reshape(fil_x,M,N)';
% y_dir=fil_x_re(:)';
% fil_xy=medfilt1(y_dir,filter_size);
% img=reshape(fil_xy,N,M)';
% figure();imshow(img);
% 
% [M,N,channel]=size(img);
% [Gmag, Gdir]=imgradient(img);
% Gmag_s=imresize(Gmag,[42,51]);
% 
% Gdir_s=imresize(Gdir,[42,51]);
% 
% nodes(:,1)=Gmag_s(:)';
% nodes(:,2)=Gdir_s(:)';
% % 间隔采样
% 
% 
% Z=linkage(nodes,'ward');
% figure();dendrogram(Z);
% c=cluster(Z,'cutoff',1.2);
% fprintf("size of cluster: %d\n",max(c));
% % figure()
% % gscatter(nodes(:,2),nodes(:,1),c);
% 
% % Define a colormap for visualization
% colormap = rand(500, 3); % Adjust the size (10) based on the maximum number of clusters
% 
% % Loop through each node in the grid
% for row = 1:42
%     for col = 1:51
%         % Get the cluster assignment for the current node
%         currentNodeCluster = c((row - 1) * 10 + col);
%         
%         % Get the color for the current cluster
%         currentColor = colormap(currentNodeCluster, :);
%         
%         % Assign the color to the pixels of the current node in the output image
%         V_merge((row - 1) * 10 + 1 : row * 10, (col - 1) * 10 + 1 : col * 10, :) = repmat(reshape(currentColor, [1, 1, 3]), [10,10,1]);
%     end
% end
% 
% % Display the resulting image
% figure();imshow(V_merge);

%% Function projectPointCloud
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
                id = floor(v) * width + floor(u);
                cloud_array(id, 1) = X(r, c);
                cloud_array(id, 2) = Y(r, c);
                cloud_array(id, 3) = z;
            end
            it = it + 1;
        end
    end
end
