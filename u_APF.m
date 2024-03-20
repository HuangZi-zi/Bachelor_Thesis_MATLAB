function [out,dir]=u_APF(img,edges)
% edges: [y，xl，xr]
[h,w,channel]=size(img);
n=2*size(edges,1);%障碍个数
%初始化车的参数
Xo=[floor(w/2),h];%起点位置
k=15;%计算引力需要的增益系数
m=3;%计算斥力的增益系数，都是自己设定的。
Po=1.5;%障碍影响距离，当障碍和车的距离大于这个距离时，斥力为0，即不受该障碍的影响。也是自己设定。

l=0.2;%步长
J=1000;%循环迭代次数
%如果不能实现预期目标，可能也与初始的增益系数，Po设置的不合适有关。
%end 

%给出障碍和目标信息，对于图像信息，先x再y
Xsum=[(edges(1,2)+edges(1,3))/2,edges(1,1);edges(:,2), edges(:,1);edges(:,3), edges(:,1)];
%这个向量是(n+1)*2维，其中第一个点[4 4]是目标位置，剩下的都是障碍的位置。
XXX=Xo;%j=1循环初始，将车的起始坐标赋给XXX
%***************初始化结束，开始主体循环******************
for j=1:J %循环开始
    goal(j,1)=XXX(1); %Goal是保存车走过的每个点的坐标。刚开始先将起点放进该向量。
    goal(j,2)=XXX(2);
    for i=1:n+1   %计算物体和障碍物、目标点的向量
         deltaX(i)=Xsum(i,1)-XXX(1);
         deltaY(i)=Xsum(i,2)-XXX(2);
         r(i)=sqrt(deltaX(i)^2+deltaY(i)^2);
    end
    Rgoal=sqrt((XXX(1)-Xsum(1,1))^2+(XXX(2)-Xsum(1,2))^2);   %路径点和目标的距离
    %目标点对路径点的引力
    Fatx=k*Rgoal*(deltaX(1)/Rgoal);
    Faty=k*Rgoal*(deltaY(1)/Rgoal);
    %各个障碍物对路径点的斥力
    for i=1:n
        if r(i+1)>Po
            Frex(i)=0;
            Frey(i)=0;
        else
            Frex(i)=-m*(1/r(i+1)-1/Po)/r(i+1)/r(i+1)*(deltaX(i+1)/r(i+1));
            Frey(i)=-m*(1/r(i+1)-1/Po)/r(i+1)/r(i+1)*(deltaY(i+1)/r(i+1));
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
        k=j;   %迭代次数
        break;
    end
end

K=j;
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

x=Xsum(2:n+1,1);
y=Xsum(2:n+1,2);

out=insertMarker(img,[x,y],'o','Color','green');
out=insertMarker(out,[X,Y],'x-mark','Color','red');
% plot(x,y,'o',X,Y,'.r');
dir=(X(5)-floor(w/2))/(Y(5)-h);
end