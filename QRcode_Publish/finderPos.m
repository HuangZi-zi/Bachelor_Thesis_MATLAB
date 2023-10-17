function [idxC, idxX, idxY] = finderPos(markers)

%% 2直線の成す角度が90度に近い点を探索
angle = zeros(3,1);
for i = 1:3
    tmp = markers;
    ref = markers(i,:);
    tmp(i,:) = [];
    % 線を表す式 : y = [x 1] * [a; b] のabを求める為に、polyfitを利用
    ab2 = polyfit([ref(1,2); tmp(1,2)], [ref(1,1); tmp(1,1)], 1);
    ab1 = polyfit([ref(1,2); tmp(2,2)], [ref(1,1); tmp(2,1)], 1);
    % ドット積を利用して角度を求める
    vect1 = [1 ab1(1)];
    vect2 = [1 ab2(1)];
    dp = dot(vect1, vect2);
    % ベクタ長の計算
    length1 = sqrt(sum(vect1.^2));
    length2 = sqrt(sum(vect2.^2));
    % 角度を算出(larger angle)
    angle(i) = 180-acos(dp/(length1*length2))*180/pi;
end
[~, idxC] = min(abs(abs(angle) - 90));

%% X軸との角度を求める
tblidx = [1 2 3];
tmp = markers;
tmp(idxC,:) = [];
tblidx(idxC) = [];
ref = markers(idxC,:);
angle = zeros(2, 1);
for i = 1:2
    % 線を表す式 : y = [x 1] * [a; b] のabを求める為に、polyfitを利用
    ab2 = polyfit([ref(1,2); tmp(i,2)], [ref(1,1); tmp(i,1)], 1);
    %ab1 = polyfit([ref(1,2); ref(1,2)], [ref(1,1); ref(1,1) + 15], 1);
    ab1(1) = 10^8; %X軸に水平な直線なので、係数は10^8程度に設定しておく
    
    % ドット積を利用して角度を求める
    vect1 = [1 ab1(1)];
    vect2 = [1 ab2(1)];
    dp = dot(vect1, vect2);
    % ベクタ長の計算
    length1 = sqrt(sum(vect1.^2));
    length2 = sqrt(sum(vect2.^2));
    % 角度を算出(larger angle)
    angle(i) = acos(dp/(length1*length2))*180/pi;
end
idx = angle > 90;
angle(idx) = 180 - angle(idx);

[~, i] = min(abs(angle));
idxX = tblidx(i);
tblidx(i) = [];
idxY = tblidx;

end

%%
% Copyright 2019 The MathWorks, Inc.