img=imread("Resource\depth_simu.png");
img=double(im2gray(img));
[G,E]=initgraph(img);

%% 



%%
function [output1,output2] = initgraph(img)
[M,N,channel]=size(img);
H=16;
W=16;

numNodesRows = floor(M / H);
numNodesCols = floor(N / W);

V=cell(numNodesRows, numNodesCols);% nodes

E=[];% edge
for i=1:numNodesRows
    for j=1:numNodesCols
        v=img((i-1)*H+1:min(i*H,M),(j-1)*W+1:min(j*W,N));
        if rejectnode(v)
            v=0;%zeros(H,W);
        end
        V{i,j}=v;
    end
end

for i=2:M/H
    for j=2:N/M
        if rejectedge(V{i,j-1},V{i,j},V{i,j+1})
            E=union(E,[V(i,j-1) V(i,j);V(i,j) V(i,j+1)]);
        elseif rejectedge(V(i-1,j),V(i,j),V(i+1,j))
            E=union(E,[V(i-1,j) V(i,j);V(i,j) V(i+1,j)]);
        end
    end
end
output1=V;
output2=E;
end

%%
function output=rejectnode(v)
[H,W,channel]=size(v);
gmag=imgradient(v);
Tmse=0.1;
count_miss=0;
for i=1:H
    for j=1:W
        if v(i,j)==0 % missing data
            count_miss=count_miss+1;
        elseif gmag(i,j)> 20 % depth-discontinuous
            output=true;
            return;
        end
    end
end
if count_miss>250
    output=true;
    return;
elseif MSE(v)>Tmse
    output=true;
    return;
else
    output=false;
end
end

%%
function output=rejectedge(va,vb,vc)
Tang=10;
normal_va=fitPlane_SVD(va);
normal_vc=fitPlane_SVD(vc);
sigma = acos(dot(normal_va,normal_vc)/(norm(normal_va)*norm(normal_vc)));%弧度制
normal_diff=abs(sigma/pi*180); %换算成角度

if va==0||vb==0||vc==0
    output=true;
    return;
elseif normal_diff>Tang
    output=true;
    return;
else
    output=false;
end
end

function [normal, MSE] = fitPlane_SVD(v)
% 功能：利用SVD拟合平面
% 输入：data   - 原始数据(m*3)    
% 输出：planes - 拟合所得平面参数 
[w,h,channel]=size(v);
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

x = 1:w;
y = 1:h;
[X,Y]= meshgrid (x,y);
ZFit = -(planes(4) + planes(1) * X + planes(2) * Y)/planes(3);
error=v-ZFit;
MSE=mse(error(:));
end

%%
function outcome=MSE(v)
if v==0
    outcome=inf;
    return;
else 
    [~,outcome]=fitPlane_SVD(v);
end
end

%%
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
        if |v|>=Tnum
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
            E=union(E,)

    
end
end

%%





%% 

