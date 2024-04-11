function [out,dir_this_time,dir_next_time,weigh]=u_APF(img,edges,v)
% edges: [y，xl，xr]
[h,w,~]=size(img);
nodesize=edges(2,1)-edges(1,1);
n=2*size(edges,1);%障碍个数
%初始化车的参数
Xo=[floor(w/2),h];%起点位置
k=2;%计算引力需要的增益系数
m=400;%计算斥力的增益系数
Po=100;%障碍影响距离，当障碍和车的距离大于这个距离时，斥力为0，即不受该障碍的影响

l=0.5;%步长
J=1000;%循环迭代次数

%给出障碍和目标信息，对于图像信息，先x再y
destination=(edges(4,2)+edges(4,3)+edges(5,2)+edges(5,3)+edges(6,2)+edges(6,3))/6;% 三行节点求平均
Xsum=[destination,edges(5,1);edges(:,2), edges(:,1);edges(:,3), edges(:,1)];
%这个向量是(n+1)*2维，其中第一个点是目标位置，剩下的都是障碍的位置。
XXX=Xo;%XXX表示车的当前位置

goal=zeros(J,2);
deltaX=zeros(n+1,1);
deltaY=zeros(n+1,1);
r=zeros(n+1,1);
Frex=zeros(n+1,1);
Frey=zeros(n+1,1);
%***************初始化结束，开始主体循环******************
for j=1:J %循环开始
    goal(j,1)=XXX(1); %Goal是保存车走过的每个点的坐标。刚开始先将起点放进该向量。
    goal(j,2)=XXX(2);
    for i=1:n+1   %计算物体和障碍物、目标点的向量
         deltaX(i)=Xsum(i,1)-XXX(1);
         deltaY(i)=Xsum(i,2)-XXX(2);
         r(i)=sqrt(deltaX(i)^2+deltaY(i)^2);
    end
    %Rgoal=sqrt((XXX(1)-Xsum(1,1))^2+(XXX(2)-Xsum(1,2))^2);   %路径点和目标的距离
    %目标点对路径点的引力
    Fatx=k*deltaX(1);
    Faty=k*deltaY(1);
    %各个障碍物对路径点的斥力
    for i=1:n
        if r(i+1)>Po
            Frex(i)=0;
            Frey(i)=0;
        else
            Frex(i)=-m*deltaX(i+1)/(r(i+1)*r(i+1));
            Frey(i)=-m*deltaY(i+1)/(r(i+1)*r(i+1));
        end
    end
    %计算合力
    Fsumx=Fatx+sum(Frex);
    Fsumy=Faty+sum(Frey);
    F=sqrt(Fsumx^2+Fsumy^2);
    %求解下一个路径点
    Xnext(1)=(XXX(1)+l*Fsumx/F);   %式子中的l是步长
    Xnext(2)=(XXX(2)+l*Fsumy/F);
    XXX=Xnext;
    
    if (sqrt((XXX(1)-Xsum(1,1))^2+(XXX(2)-Xsum(1,2))^2)<0.1)   %当物体接近目标点时
%         K=j;   %迭代次数
        break;
    end
end
K=j;

weigh=2^(-sqrt((XXX(1)-Xsum(1,1))^2+(XXX(2)-Xsum(1,2))^2));

goal(K,1)=Xsum(1,1);%把路径向量的最后一个点赋值为目标
goal(K,2)=Xsum(1,2);

%***********************************画出障碍，起点，目标，路径点*************************
% figure(3);imshow(img);
% hold on;
%画出路径
X=goal(:,1);
Y=goal(:,2);
%路径向量Goal是二维数组,X,Y分别是数组的x,y元素的集合，是两个一维数组。
%x=[1 3 4 3 6 5.5 8];%障碍的x坐标
%y=[1.2 2.5 4.5 6 2 5.5 8.5];

% 画出障碍物
x=Xsum(2:n+1,1);
y=Xsum(2:n+1,2);

out=insertMarker(img,[x,y],'o','Color','green');
out=insertMarker(out,[X,Y],'x-mark','Color','red');
% plot(x,y,'o',X,Y,'.r');

line1=abs(Y-(h-nodesize));
[~, index1] = min(line1(:));
% dir_this_time=2*v./(1+exp(-(1/v*10)*(X(index1)-w/2)))-v;
dir_this_time=200./(1+exp(-(1/50)*(X(index1)-w/2)))-100;
% disp(X(index1));
line2=abs(Y-(h-2*nodesize));
[~, index2] = min(line2(:));

dir_next_time=200./(1+exp(-(1/50)*(X(index2)-X(index1))))-100;

end