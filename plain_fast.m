clear;
img=imread("Resource\depth_simu.png");
img=double(im2gray(img));
% [G,E]=initgraph(img);
[M,N,channel]=size(img);
H=16;
W=16;

Tmse=100;


% data=zeros(M*N,3);
% for i=1:M
%     for j=1:N
%         data((i-1)*N+j,:)=[i;j;img(i,j)];
%     end
% end

numNodesRows = floor(M / H);
numNodesCols = floor(N / W);

V_pixel=cell(numNodesRows, numNodesCols);% nodes
V_normal=cell(numNodesRows, numNodesCols);
V_angle=zeros(numNodesRows, numNodesCols);
V_plane=cell(numNodesRows, numNodesCols);
V_MSE=zeros(numNodesRows, numNodesCols);

V_merge=zeros(M,N);

merged=cell(1,numNodesRows*numNodesCols/2);

z_mean=zeros(numNodesRows, numNodesCols);

% 每个node中的相对座标
points_r(:,1)=repelem(1:H,W)';%x
points_r(:,2)=repmat(1:W,1,W)';%y
% 归化为均值为0的序列后每个node中的座标
points_c(:,1)=points_r(:,1)-(1+W)/2;
points_c(:,2)=points_r(:,2)-(1+H)/2;

% v=struct('pixels', [], 'normal', [], 'plane', [],'v_MSE',[]);
E=zeros(H*W,4);% edge

for i=1:numNodesRows
    for j=1:numNodesCols
        V_pixel{i,j}=img((i-1)*H+1:min(i*H,M),(j-1)*W+1:min(j*W,N));
       
        points(:,1)=points_r(:,1)+H*(i-1);%x
        points(:,2)=points_r(:,2)+W*(j-1);%y
        points(:,3)=(V_pixel{i,j}(:));%z

        z_mean(i,j)=mean(points(:,3));

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
        V_normal{i,j}=eigenvector(index(1),:)./norm(eigenvector(index(1),:));
        gmag=imgradient(V_pixel{i,j});
        
        if  any([z_mean(i,j)<2,max(max(gmag))>200,sorted_enigenvalue(1)>Tmse]) %missing %discontinue %large mse
            V_MSE(i,j)=inf;
            V_pixel{i,j}=[];
            E=union(E,[i,j-1,i,j ; i,j,i,j+1 ; i-1,j,i,j ; i,j,i+1,j],'rows');
        else
            V_MSE(i,j)=sorted_enigenvalue(1);
        end  
    end
end

V_MSE_T=V_MSE';
heap=(V_MSE_T(:))';
[sorted_heap, index]=sort(heap);

merge_cluster=1;

for i=1:size(heap,2)
    current=index(1);
    current_i=floor((current-1)/numNodesCols)+1;
    current_j=current-(current_i-1)*numNodesCols;%rem(current,numNodesCols);
    if current_i==1||current_i==numNodesRows||current_j==1||current_j==numNodesCols
        continue
    end

    if ~isempty(V_pixel{current_i,current_j})
        points_current(:,1)=points_r(:,1)+H*(current_i-1);%x
        points_current(:,2)=points_r(:,2)+W*(current_j-1);%y
        points_current(:,3)=(V_pixel{current_i,current_j}(:));%z
%         index(1)=[];
        if isempty(merged{1,1})
            merged{1,1}=points_current; 
        else
            for j=1:merge_cluster
                older=union(merged{1,j},points_current,"rows");
                [outcome,newMse]=trymerge(older,Tmse);
                if outcome
                    merge_cluster=j;
                    V_merge(points_current(:,1),points_current(:,2))=merge_cluster;
                    merged{1,merge_cluster}=older;
                    index=setdiff(index,current,'stable');
                    break;
                elseif j==merge_cluster
                    % 如果没有break，说明不能merge进之前的块里
                    merge_cluster=merge_cluster+1;
                    merged{1,merge_cluster}=points_current;
                end
            end

        end
    else
        continue
    end

    
    if ~isempty(V_pixel{current_i,current_j-1}) && ismember(current-1,index)
        points_left(:,1)=points_r(:,1)+H*(current_i-1);%x
        points_left(:,2)=points_r(:,2)+W*(current_j-1-1);%y
        points_left(:,3)=(V_pixel{current_i,current_j-1}(:));%z
        merged_left=union(merged{1,merge_cluster},points_left,"rows");
        [outcome,newMse]=trymerge(merged_left,Tmse);
        if outcome
            V_merge(points_left(:,1),points_left(:,2))=merge_cluster;
            merged{1,merge_cluster}=merged_left;
            index=setdiff(index,current-1,'stable');
        else
            merge_cluster=merge_cluster+1;
            merged{1,merge_cluster}=points_current;
        end
    end
    
    if ~isempty(V_pixel{current_i,current_j+1}) && ismember(current+1,index)
        points_right(:,1)=points_r(:,1)+H*(current_i-1);%x
        points_right(:,2)=points_r(:,2)+W*(current_j-1+1);%y
        points_right(:,3)=(V_pixel{current_i,current_j+1}(:));%z
        merged_right=union(merged{1,merge_cluster},points_right,"rows");
        [outcome,newMse]=trymerge(merged_right,Tmse);
        if outcome
            V_merge(points_right(:,1),points_right(:,2))=merge_cluster;
            merged{1,merge_cluster}=merged_right;
            index=setdiff(index,current+1,'stable');
        else
            merge_cluster=merge_cluster+1;
            merged{1,merge_cluster}=points_current;
        end
    end

    if ~isempty(V_pixel{current_i-1,current_j}) && ismember(current-numNodesCols,index)
        points_up(:,1)=points_r(:,1)+H*(current_i-1-1);%x
        points_up(:,2)=points_r(:,2)+W*(current_j-1);%y
        points_up(:,3)=(V_pixel{current_i-1,current_j}(:));%z
        merged_up=union(merged{1,merge_cluster},points_up,"rows");
        [outcome,newMse]=trymerge(merged_up,Tmse);
        if outcome
            V_merge(points_up(:,1),points_up(:,2))=merge_cluster;
            merged{1,merge_cluster}=merged_up;
            index=setdiff(index,current-numNodesCols,'stable');
        else
            merge_cluster=merge_cluster+1;
            merged{1,merge_cluster}=points_current;
        end
    end

    if ~isempty(V_pixel{current_i+1,current_j}) && ismember(current+numNodesCols,index)
        points_down(:,1)=points_r(:,1)+H*(current_i-1+1);%x
        points_down(:,2)=points_r(:,2)+W*(current_j-1);%y
        points_down(:,3)=(V_pixel{current_i+1,current_j}(:));%z
        merged_down=union(merged{1,merge_cluster},points_down,"rows");
        [outcome,newMse]=trymerge(merged_down,Tmse);
        if outcome
            V_merge(points_down(:,1),points_down(:,2))=merge_cluster;
            merged{1,merge_cluster}=merged_down;
            index=setdiff(index,current+numNodesCols,'stable');
        else
            merge_cluster=merge_cluster+1;
            merged{1,merge_cluster}=points_current;
        end
    end
    imshow(uint8(V_merge.*100));
%     right=current+1;
%     left=current-1;
%     up=current-numNodesCols;
%     down=current+numNodesCols;  
end

%imshow(uint8(V_merge.*10));
    
function [outcome,newMSE]=trymerge(points,Tmse)
W=16;H=16;
x=points(:,1);
y=points(:,2);
z=points(:,3);

pc(:,1)=x-mean(x);
pc(:,2)=y-mean(y);
pc(:,3)=z-mean(z);

sigma=pc'*pc./(W*H);
        [eigenvector,eigenvalue]=eig(sigma);
        eigenvalue=diag(eigenvalue);
        [sorted_enigenvalue, index]=sort(eigenvalue);
        if sorted_enigenvalue(1)<Tmse
            outcome=true;
            newMSE=sorted_enigenvalue(1);
        else
            outcome=false;
            newMSE=inf;
        end

end



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
function output=rejectedge(va,vb,vc,na,nc)
% 123输入像素 45输入法向量
Tang=10;
normal_va=na;
normal_vc=nc;
sigma = acos(dot(normal_va,normal_vc)/(norm(normal_va)*norm(normal_vc)));%弧度制
normal_diff=abs(sigma/pi*180); %换算成角度

% if va==0||vb==0||vc==0
if isempty(va)||isempty(vb)||isempty(vc)
    output=true;
    return;
elseif normal_diff>Tang
    output=true;
    return;
else
    output=false;
end
end

%% 
function [normal, planes] = fitPlane_SVD(v)
% 功能：利用SVD拟合平面
% 输入：data   - 原始数据(m*3)    
% 输出：planes - 拟合所得平面参数 
% v=input.pixels
[w,h,channel]=size(v);
data=zeros(w*h,3);
for i=1:w
    for j=1:h
        data((i-1)*h+j,:)=[i;j;v(i,j)];
    end
end
points = data(:,1:3);
% 去质心
M = points - repmat(mean(points),size(points,1),1);
% 奇异值分解
[~,~,V] = svd(M,0);
% 最小特征值对应的特征向量为法向量
normal = V(:,3)';
% 平面参数标准化
dtmp = mean(points*normal');
planes(1:3) = normal'*sign(dtmp);
planes(4) = -dtmp*sign(dtmp);

% x = 1:w;
% y = 1:h;
% [X,Y]= meshgrid (x,y);
% ZFit = -(planes(4) + planes(1) * X + planes(2) * Y)/planes(3);
% error=v-ZFit;
% v_MSE=mse(error(:));
end

% %%
% % function outcome=MSE(input)
% % v=input.pixels;
% % if v==0
% %     outcome=inf;
% %     return;
% % else 
% %     [~,outcome]=fitPlane_SVD(input);
% % end
% % end
% 
% %%
function ourput=AHCluster(V,E)
Tmse=0.1;
Tnum=160;
Q=buildminMSEheap(V);% build a min-heap data structure
B=[];
PI=[];% 平面的方程
while ~isempty(Q)
    v=popmin(Q);
    if v notin V
        continue
    end
    u_best=[];
    u_merge=[];
    for i=1:size(E,2)
        u_test=union(E(i,1),E(i,2));
        if MSE(u_test)<MSE(u_merge)
            u_best=E(i,1);
            u_merge=u_test;
        end
    end
    if MSE(u_merge)>=Tmse
        if abs(v)>=Tnum
            B=union(B,v);
            norm=fitPlane_SVD(v);
            dtmp = mean(points*norm');
            planes(1:3) = norm'*sign(dtmp);
            planes(4) = -dtmp*sign(dtmp);
            PI=union(PI,planes);
            %从E中删除E（v）
            E=setdiff(E,E(v));
            %从V中删除v
            V=setdiff(V,v);
        else
            insert(Q,u_merge);
            %E=union(E,)
        end
    end

    
end
end

%%





%% 

