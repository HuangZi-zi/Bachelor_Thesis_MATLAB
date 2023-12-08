
%%
function output = initgraph(img)
[M,N,channel]=size(img);
H=16;
W=16;

V=[];% nodes 三维点的集合
E=[];% edge
for i=1:M/H
    for j=1:N/W
        v(i,j)=img((i-1)*H:min(i*H,M),(j-1)*W:min(j*W,N));
        if rejectnode(v(i,j))
            v(i,j)=[];
        end
        V=union(V,v(i,j));
    end
end

for i=1:M/H
    for j=1:N/M
        if rejectedge(V(i,j-1),V(i,j),V(i,j+1))
            E=union(E,[V(i,j-1) V(i,j);V(i,j) V(i,j+1)]);
        elseif rejectedge(V(i-1,j),V(i,j),V(i+1,j))
            E=union(E,[V(i-1,j) V(i,j);V(i,j) V(i+1,j)]);
        end
    end
end
output=[V;E];
end

%%
function output=rejectnode(v)
[H,W,channel]=size(v);
gmag=imgradient(v);
Tmse=0.1;
for i=1:H
    for j=1:W
        if isempty(v(i,j)) % missing data
            output=true;
            return;
        elseif gmag(i,j)> 20 % depth-discontinuous
            output=true;
            return;
        end
    end
end
if MSE(v)>Tmse
    output=true;
    return;
else
    output=false;
end
end

%%
function output=rejectedge(va,vb,vc)
Tang=5;
normal_va=fitPlane_SVD(va);
normal_vc=fitPlane_SVD(vc);
sigma = acos(dot(normal_va,normal_vc)/(norm(normal_va)*norm(normal_vc)));%弧度制
normal_diff=abs(sigma/pi*180); %换算成角度

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

function planes = fitPlane_SVD(data)
% 功能：利用SVD拟合平面
% 输入：data   - 原始数据(m*3)    
% 输出：planes - 拟合所得平面参数 
points = data(:,1:3);
% 去质心
M = points - repmat(mean(points),size(points,1),1);
% 奇异值分解
[~,~,V] = svd(M,0);
% 最小特征值对应的特征向量为法向量
normal = V(:,3)';
% 平面参数标准化
% dtmp = mean(points*normal');
% planes(1:3) = normal'*sign(dtmp);
% planes(4) = -dtmp*sign(dtmp);
planes=normal;
end

