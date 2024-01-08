function [outputArg1,outputArg2] = u_plane_ahc(pointx,pointy,pointz,index,width,height)

N=width;
M=height;

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
V_normal=zeros(numNodesRows, numNodesCols, 3);
V_angle=zeros(numNodesRows, numNodesCols);
V_plane=cell(numNodesRows, numNodesCols);
V_MSE=zeros(numNodesRows, numNodesCols);

V_merge=zeros(M,N,3);

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

% 主成分分析
for i=1:numNodesRows
    for j=1:numNodesCols
        V_pixel{i,j}=img((i-1)*H+1:min(i*H,M),(j-1)*W+1:min(j*W,N));
       
        points(:,1)=points_r(:,1)+H*(i-1);%x
        points(:,2)=points_r(:,2)+W*(j-1);%y
        points(:,3)=(V_pixel{i,j}(:));%z

        z_mean(i,j)=mean(points(:,3));

        points_c(:,3)=points(:,3)-z_mean(i,j);

        sigma=points_c'*points_c./(W*H);

        [eigenvector,eigenvalue]=eig(sigma);
        eigenvalue=diag(eigenvalue);
        [sorted_enigenvalue, index]=sort(eigenvalue);
        normal_tempt=eigenvector(index(1),:)./norm(eigenvector(index(1),:));
        V_normal(i,j,:)=reshape(normal_tempt, [1, 1, 3]);
        V_angle(i,j)=acos(dot(normal_tempt,[1,1,1])/(norm(normal_tempt)*norm([1,1,1])))/pi*180;
        gmag=imgradient(V_pixel{i,j});
        
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
c=cluster(Z,'cutoff',80,'criterion','distance');
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
    


end